#!/usr/bin/python

# Generic Methods
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_lrbridge
Description:  Long read bridging methods library
Version:      0.0
Last Edit:    19/04/21
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains methods, modified from Diploidocus, to underpin GapSpanner, LocusFixer and SynBad reassembly
    of genome assembly regions. In each case, an RJE Object will need to be given to the method as self. This will enable
    easy incorporation of existing methods, with a quick update from self.method() calls to rje_lrbridge.method(self)
    calls.

Commandline:
    This module is not for standalone running and has no commandline options (including 'help'). All options are handled
    by the parent modules.

Uses general modules: copy, os, string, sys
Uses RJE modules: rje
Other modules needed: rje_py2, rje_py3
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
### <1> ### Region span and reassembly methods                                                                          #
#########################################################################################################################
#i# The gapSpan method has been directly copied from Diploidocus
def gapSpan(self): ### Performs assembly gap read check and spanning analysis.
    '''
    ### ~ Assembly gap read-spanning analysis [runmode=gapspan] ~ ###

    This mode first identifies all the gaps in an assembly (`seqin=FILE`) (using SeqList `gapstats` or `$SEQBASE.gaps.tdt` if pre-
    existing) and then runs the read spanning analysis (`runmode=regcheck`) with `regcnv=F`. Long read data, given
    with the `reads=FILELIST` and `readtype=LIST` options, are mapped onto the assembly using minimap2 to generate a PAF file.
    This is then parsed and reads spanning each gap are identified based on their positions and the target start and end positions in the PAF file.
    In addition to absolute spanning of regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
    sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
    using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
    due to sequence length constraints.

    Spanning `spanid` output is also generated for each gap and saved in `$BASEFILE_spanid`. Each gap will be named:
    `seqname.start-end`.
    '''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        db = self.db()
        mainbase = self.baseFile()
        basefile = self.baseFile(strip_path=True)
        seqbase =  rje.baseFile(self.getStr('SeqIn'),strip_path=True)
        ## ~ [1a] ~ Check or create gapstats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=seqbase,log=self.log):
            self.cmd_list.append('gapstats')    # Try to run automatically if possible
            seqin = self.seqinObj()
            if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=seqbase,log=self.log):
                seqin.setBool({'Raw':False,'GapStats':True,'DNA':True})
                seqin.str['SeqType'] = 'dna'
                seqin.summarise()
        if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=seqbase,log=None):
            raise ValueError('Problem finding/generating %s.gaps.tdt' % seqbase)
        gapfile = '%s.gaps.tdt' % seqbase
        checkfile = '%s.checkpos.tdt' % basefile
        ## ~ [1b] ~ Setup PAF CheckFlanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not self.list['CheckFlanks']:
            self.list['CheckFlanks'] = [0]
            self.cmd_list.append('checkflanks=0')
        minflank = self.list['CheckFlanks'][0]
        minspan = 'Span%d' % minflank
        maxflank = self.list['CheckFlanks'][-1]

        ### ~ [2] ~ Load Gaps and run checkpos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cdb = None
        if rje.checkForFiles(filelist=['.checkpos.tdt'],basename=basefile,log=self.log) and not self.force():
            cdb = db.addTable(checkfile,mainkeys=['seqname','start','end'],name='checkpos',ignore=[],expect=True)
            cdb.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            if minspan not in cdb.fields() and 'Span0' in cdb.fields():
                self.printLog('#SPAN','Cannot find "{0}" field in {1}: will use Span0'.format(minspan,checkfile))
                minspan = 'Span0'
            elif minspan not in cdb.fields():
                self.printLog('#SPAN','Cannot find "{0}" field in {1}: will regenerate'.format(minspan,checkfile))
                cdb = None
        if not cdb:
            cdb = db.addTable(gapfile,mainkeys=['seqname','start','end'],name='gap_pos',ignore=[],expect=True)
            cdb.dataFormat({'seqlen':'int','start':'int','end':'int'})
            cdb.makeField(formula='#seqname#.#start#-#end#',fieldname='gap')
            cdb.saveToFile()
            self.setStr({'RegCheck':cdb.saveToFileName()})
            self.list['CheckFields'] = ['seqname','start','end']
            checkcmd = ['checkfields=seqname,start,end','spanid=gap']
            ## ~ [2a] ~ Run Checkpos and SpanID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            vfile = self.getStr('RegCheck')
            pfile = self.getPAFFile()   #baseFile() + 'checkpos.paf'
            pafcmd = self.cmd_list + ['checkpos={}'.format(vfile),'pafin={}'.format(pfile)] + checkcmd
            paf = rje_paf.PAF(self.log, pafcmd)
            cdb = paf.checkPos(save=False)
            cdb.dataFormat({'MaxFlank3':'int','seqlen':'int','end':'int'})
            for centry in cdb.entries(): centry['MaxFlank3'] = centry['seqlen'] - centry['end']
            cdb.rename('checkpos')
            cdb.saveToFile(backup=False)
            #!# Clean up gap_pos file #!#

        ### ~ [3] ~ Reassemble gaps (runmode=gapass) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        asslist = []    # List of assembly commands to fork out and outputs (acmd, output, assembly)
        reassemble = self.getStrLC('RunMode') in ['gapass','gapfill']
        assembler = 'flye' # self.getStr('Assembler')
        filldir = '%s_gapfill/' % basefile
        assfile = '{0}{1}.assembledgaps.fasta'.format(filldir,basefile)
        if self.getStrLC('RunMode') in ['gapfill'] and rje.exists(assfile) and not self.force():
            self.printLog('#GAPASS','{0} found (force=F): skipping assembly'.format(assfile))
            reassemble = False
        if reassemble:
            try:
                vcheck = rje.split(os.popen('{0} --version'.format(assembler)).read())[0]
                self.printLog('#PROG','{0} version: {1}'.format(assembler,vcheck))
            except:
                raise ValueError('Assembler check error - failed to run: {0} --version'.format(assembler))
        if reassemble:
            rtype = self.list['ReadType'][0]
            if rtype in ['pacbio','pac']: rtype = 'pb'
            if rtype in ['hifi','ccs']: rtype = 'hifi'
            if rtype not in ['ont','pb','hifi']:
                self.warnLog('Read Type "%s" not recognised (pb/ont/hifi): check readtype=LIST. Will use "ont".' % rtype)
                rtype = 'ont'
            if assembler == 'flye':
                if rtype == 'ont': rtype = 'nano-raw'
                elif rtype == 'hifi': rtype = 'pacbio-hifi'
                else: rtype = 'pacbio-raw'
            bamfile = self.getBamFile()
            iddir = '%s_spanid/' % basefile
            assdir = '%s_gapassemble/' % basefile
            rje.mkDir(self,assdir)
            nonx = 0; assx = 0; failx = 0; skipx = 0
            for centry in cdb.entries():
                spanner = centry['gap']
                target = '{0}{1}.{2}.assembly.fasta'.format(assdir,basefile,spanner)
                if rje.exists(target) and not self.force():
                    self.printLog('#GAPASS','{0} found (force=F): skipped.'.format(target))
                    skipx = 0
                    continue
                spout = '%s%s.%s.span.id' % (iddir,basefile,spanner)
                ## ~ [3a] ~ Check for read ID files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not rje.exists(spout):
                    if centry[minspan] > 0:
                        self.warnLog('No {0} file but {1}={2}. May need to delete {3} and regenerate.'.format(spout,minspan,centry[minspan],checkfile))
                    elif self.v() > 1: self.printLog('No %s file' % spout)
                    nonx += 1
                    continue
                ## ~ [3b] ~ Extract the reads from the BAM file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #!# Convert this to a forked process too?
                alog = '{0}{1}.{2}.assemble.log'.format(assdir,basefile,spanner)
                self.printLog('#BAM','Extracting gap-spanning reads from BAM file.')
                # Exclude secondary and umapped reads
                region = '{0}:{1}-{2}'.format(centry['seqname'],max(centry['start']-maxflank,1),centry['end']+maxflank)
                idbam = '%s%s.%s.bam' % (assdir,basefile,spanner)
                fasout = '%s%s.%s.reads.fasta' % (assdir,basefile,spanner)
                if self.force() or not rje.exists(fasout):
                    gap2bam = 'samtools view -h -F 0x104 {0} {1} | grep -e @ -f {2} | samtools view -h -b -o {3} -'.format(bamfile,region,spout,idbam)
                    logline = self.loggedSysCall(gap2bam,alog,append=True,stderr=False,nologline='No stdout from samtools view',threaded=False)
                    bam2fas = 'samtools fasta {0} > {1}'.format(idbam,fasout)
                    logline = self.loggedSysCall(bam2fas,alog,append=True,stderr=False,nologline='No stdout from samtools fasta',threaded=False)
                gensize = maxflank*2
                if rje.exists(fasout) and open(fasout,'r').read():
                    fascmd = ['dna=T','raw=T','gapstat=F','summarise=F','autoload=T','seqmode=list','seqin={0}'.format(fasout)]
                    fasta = rje_seqlist.SeqList(self.log,self.cmd_list+fascmd)
                    if fasta.seqNum() < 1:
                        self.warnLog('{0} has no reads despite {1} spanning read ID output'.format(fasout,spout))
                        continue
                    fasdat = fasta.summarise(sumdb=False,save=False)
                    gensize = int(fasdat['MeanLength'])
                else:
                    self.warnLog('{0} has no reads despite {1} spanning read ID output'.format(fasout,spout))
                    continue
                ## ~ [3c] Assemble reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Adding forking
                #i# asslist = []    # List of assembly commands to fork out and outputs (acmd, output, assembly)
                if assembler == 'flye':
                    acmd = 'flye --{0} {1} --out-dir {2}{3}.{4}_flye --genome-size {5} --threads {6}'.format(rtype,fasout,assdir,basefile,spanner,gensize,self.getInt('SubForks'))
                    assembly = '{0}{1}.{2}_flye/assembly.fasta'.format(assdir,basefile,spanner)
                if self.v() > 0:
                    acmd = '{0} 2>&1 | tee -a {1}'.format(acmd,alog)
                else:
                    acmd = '{0} 2>&1 >> {1}'.format(acmd,alog)
                asslist.append((acmd,assembly,'{0}{1}.{2}.assembly.fasta'.format(assdir,basefile,spanner)))
                continue
                #!# Old code:
                assembly = 'Assembly'
                self.printLog('#GAPASS','Assembling gap-spanning reads for {0}.'.format(spanner))
                if assembler == 'flye':
                    acmd = 'flye --{0} {1} --out-dir {2}{3}.{4}_flye --genome-size {5} --threads {6}'.format(rtype,fasout,assdir,basefile,spanner,gensize,self.getInt('SubForks'))
                    logline = self.loggedSysCall(acmd,alog,append=True,nologline='No stdout from flye',threaded=False)
                    assembly = '{0}{1}.{2}_flye/assembly.fasta'.format(assdir,basefile,spanner)
                if rje.exists(assembly):
                    self.printLog('#GAPASS','{0} generated -> {1}'.format(assembly,target))
                    rje.fileTransfer(assembly,target,deletefrom=False,append=False)
                    assx += 1
                else:
                    self.printLog('#GAPASS','No {0} generated.'.format(assembly))
                    failx += 1
            #self.printLog('#GAPASS','{0} Gap assemblies generated; {1} failed; {2} gaps without read IDs'.format(rje.iStr(assx),rje.iStr(failx),rje.iStr(nonx)))
            self.printLog('#GAPASS','{0} gaps without read IDs'.format(rje.iStr(nonx)))
            if skipx:
                self.printLog('#GAPASS','{0} existing assemblies skipped (force=F).'.format(rje.iStr(skipx)))

            ### ~ [4] ~ Fork out assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#FORK','{0} gap assemblies to fork out'.format(rje.iLen(asslist)))
            ## ~ [4a] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            forker = self.obj['Forker']
            forker.list['ToFork'] = []
            for afork in asslist:
                forker.list['ToFork'].append(afork[0])
            self.debug('\n'.join(forker.list['ToFork']))
            ## ~ [4b] Fork out assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if forker.list['ToFork']:
                if self.getNum('Forks') < 1:
                    #i# Warn lack of forking
                    self.printLog('#FORK','Note: program can be accelerated using forks=INT.')
                    for acmd in forker.list['ToFork']:
                        alog = afork[2][:-6] + '.log'
                        logline = self.loggedSysCall(acmd,alog,append=True,nologline='No stdout from flye',threaded=False)
                        #self.printLog('#SYS',forkcmd)
                        #os.system(forkcmd)
                else:
                    self.printLog('#FORK','Forking assemblies using {0} x {1} threads'.format(self.getNum('Forks'),self.getNum('SubForks')))
                    if forker.run():
                        self.printLog('#FORK','Forking of assemblies completed.')
                    else:
                        try:
                            self.errorLog('Assembly forking did not complete',printerror=False,quitchoice=True)
                        except:
                            raise RuntimeError('Assembly forking did not complete')
            ## ~ [4c] Process assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            assx = 0; failx = 0
            for afork in asslist:
                assembly = afork[1]
                target = afork[2]
                if rje.exists(assembly):
                    self.printLog('#GAPASS','{0} generated -> {1}'.format(assembly,target))
                    rje.fileTransfer(assembly,target,deletefrom=False,append=False)
                    assx += 1
                else:
                    self.printLog('#FAIL','No {0} generated.'.format(assembly))
                    failx += 1
            self.printLog('#GAPASS','{0} Gap assemblies generated; {1} failed'.format(rje.iStr(assx),rje.iStr(failx)))

        ### ~ [5] ~ Fill gaps (runmode=gapfill) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #i# The idea here is to use GABLAM-style PAF mapping to map the assembled regions back onto the assembly
        if self.getStrLC('RunMode') in ['gapfill']:
            ## ~ [4a] ~ Compile the assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            assdir = '%s_gapassemble/' % basefile
            filldir = '%s_gapfill/' % basefile
            rje.mkDir(self,filldir)
            assfile = '{0}{1}.assembledgaps.fasta'.format(filldir,basefile)
            if self.force() or not rje.exists(assfile):
                rje.backup(self,assfile)
                ASSFILE = open(assfile,'w')
                nonx = 0; assx = 0; nullx = 0; ctgx = 0
                for centry in cdb.entries():
                    spanner = centry['gap']
                    assembly = '{0}{1}.{2}.assembly.fasta'.format(assdir,basefile,spanner)
                    if not rje.exists(assembly): nonx += 1; continue
                    aseq = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=list','autoload=F','seqin={0}'.format(assembly)])
                    aseq.loadSeq()
                    if aseq.list['Seq']:
                        sx = 1; assx += 1
                        for (name, sequence) in aseq.list['Seq']:
                            ASSFILE.write('>{0}.ctg{1} {2}\n{3}\n'.format(spanner,sx,name,sequence))
                            sx += 1; ctgx += 1
                    else: nullx += 1
                ASSFILE.close()
                self.printLog('#GAPCTG','Collated {0} assembled contigs for {1} gaps; {2} assemblies w/o contigs; {3} w/o assemblies'.format(rje.iStr(ctgx),rje.iStr(assx),rje.iStr(nullx),rje.iStr(nonx)))
            else:
                self.printLog('#GAPCTG','Collated assembled contigs found (force=F): {0}'.format(assfile))
            ## ~ [4b] ~ Map gap assemblies to genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fillbase = '%s%s.gapctg' % (filldir,self.baseFile())
            pafin = '%s.paf' % (fillbase)
            paf = rje_paf.PAF(self.log, self.cmd_list+['pafin=%s' % pafin,'%s%s.gapctg' % (filldir,self.baseFile()),'seqin=%s' % assfile,'reference=%s' % self.getStr('SeqIn')])
            defaults = paf_defaults
            paf.dict['MapOpt'] = rje.combineDict(defaults,paf.dict['MapOpt'],overwrite=True)
            paf.setup()
            if not rje.exists(pafin) or self.force():
                rje.backup(self,pafin)
                paf.minimap2()
                #!# Add localaln=T to PAF #!#
            #?# Should paf.db() basename be changed to pafin basefile for the initial output?
            ## ~ [4c] ~ GapFiller local hits table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gapfiller = '%s.gapfiller.tdt' % (fillbase)
            paf.basefile(fillbase)
            if rje.exists(gapfiller) and not self.force():
                udb = paf.db().addTable(gapfiller,mainkeys=['Hit','SbjStart','SbjEnd'],name='gapfiller',ignore=[],expect=True)
                udb.dataFormat({'QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','QryLen':'int','SbjLen':'int'})
            else:
                paf.run()
                paf.db('hitunique').rename('gapfiller')     # ['Qry','QryStart','QryEnd']
                paf.db('gapfiller').dropFields(['nn','cs'])
                for entry in paf.db('gapfiller').entries():
                    if entry['SbjStart'] > entry['SbjEnd']:
                        (entry['SbjStart'],entry['SbjEnd']) = (entry['SbjEnd'],entry['SbjStart'])
                paf.db('gapfiller').newKey(['Hit','SbjStart','SbjEnd'])
                paf.db('gapfiller').saveToFile()
                #  Qry     Hit     AlnNum  BitScore        Expect  Length  Identity        Positives       QryStart        QryEnd  SbjStart        SbjEnd  QryLen  SbjLen  Strand
            self.printLog('#GAPMAP','%s%s.gapctg.* output generated.' % (filldir,self.baseFile()))
            ## ~ [4d] ~ Identify gap fills using local hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Work through by qry position, looking for a single hit or pair of hits that spans the gap
            #i# -> For a single hit, replace the whole local hit region
            #i# -> For a pair of hits, replace the gap between the Sbj regions
            #i# => Build a list of (seqname, start, end, gapctg, start, end, strand)
            #i# .. Reverse sort of replace sequences in turn
            #i# => Output to *.gapfill.fasta
            #i# -> Report number of gaps filled and new stats
            paf.basefile(mainbase)
            udb = paf.db('gapfiller')
            preventry = None
            gapfills = []   # Tuples of (seqname, start, end, gapctg, start, end, strand)
            filled = []
            #i# Filter on local alignment length (and identity?)
            minloclen = self.getInt('MinLocLen')
            minlocid = self.getPerc('MinLocID')
            if minloclen > 1: udb.dropEntries('Length<%d' % minloclen)
            if minlocid > 0:
                udb.makeField('100.0*Identity/Length','LocID')
                udb.dropEntries('LocID<%s' % minlocid)
                udb.dropField('LocID')
            #i# Identify gaps to fill
            ex = 0.0; etot = udb.entryNum()
            for entry in udb.entries(sorted=True):
                self.progLog('\r#FILL','Identifying gap-filling re-assemblies: %.2f%%' % (ex/etot)); ex += 100.0
                gapid = '.'.join(entry['Qry'].split('.')[:-1])  # PTEXCHR1.01.10697073-10697172.ctg1
                if gapid in filled: continue
                #i# Check for matching sequences
                if not entry['Qry'].startswith('{0}.'.format(entry['Hit'])): continue
                #i# Check preventry
                if preventry and preventry['Qry'] != entry['Qry']: preventry = None
                #i# Check for single hit spanning gap
                #if entry['Strand'] == '-': (entry['SbjStart'],entry['SbjEnd']) = (entry['SbjEnd'],entry['SbjStart'])
                [gapx,gapy] = entry['Qry'].split('.')[-2].split('-')
                gapx = int(gapx)
                gapy = int(gapy)
                if entry['SbjStart'] < gapx and entry['SbjEnd'] > gapy:
                    filled.append(gapid)
                    gapfills.append((entry['Hit'],entry['SbjStart'],entry['SbjEnd'],entry['Qry'],entry['QryStart'],entry['QryEnd'],entry['Strand'],'span'))
                    preventry = entry
                    continue
                if not preventry:
                    preventry = entry
                    continue
                #i# Check for dual hit spanning gap
                if preventry['SbjEnd'] <= gapx and entry['SbjStart'] >= gapy and preventry['Strand'] == entry['Strand']:
                    if (entry['Strand'] == '+' and entry['SbjStart'] > preventry['SbjEnd']):
                        filled.append(gapid)
                        gapfills.append((entry['Hit'],preventry['SbjEnd'],entry['SbjStart'],entry['Qry'],preventry['QryEnd'],entry['QryStart'],entry['Strand'],'flank'))
                    elif (entry['Strand'] == '-' and preventry['SbjStart'] > entry['SbjEnd']):
                        filled.append(gapid)
                        gapfills.append((entry['Hit'],preventry['SbjEnd'],entry['SbjStart'],entry['Qry'],entry['QryEnd'],preventry['QryStart'],entry['Strand'],'flank'))
                preventry = entry
            self.progLog('\r#FILL','Identified %s gap-filling re-assemblies from %s.' % (rje.iLen(gapfills),gapfiller))
            ## ~ [4e] ~ GapFill using local hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add a check for overlapping replacements? #!#
            gapfills.sort()
            gapshift = {}
            gdb = db.addEmptyTable('gapfill',['seqname','start','end','gapctg','ctgstart','ctgend','strand','gaplen','type','newname','newstart','newend','newchunklen'],['seqname','start','end'],log=False)
            for filldat in gapfills:
                fentry = {'seqname':filldat[0],'start':filldat[1],'end':filldat[2],'gapctg':filldat[3],'ctgstart':filldat[4],'ctgend':filldat[5],'strand':filldat[6],'type':filldat[7]}
                if fentry['seqname'] not in gapshift: gapshift[fentry['seqname']] = 0
                fentry['newname'] = fentry['seqname']
                fentry['newstart'] = fentry['start'] + gapshift[fentry['seqname']]
                gapshift[fentry['seqname']] += (fentry['ctgend'] - fentry['ctgstart']) - (fentry['end'] - fentry['start'])
                fentry['newend'] = fentry['end'] + gapshift[fentry['seqname']]
                fentry['newchunklen'] = fentry['newend'] - fentry['newstart'] + 1
                [gapx,gapy] = fentry['gapctg'].split('.')[-2].split('-')
                fentry['gaplen'] = int(gapy) - int(gapx) + 1
                gdb.dict['Data'][(filldat[0],filldat[1],filldat[2])] = fentry
            gdb.index('seqname')
            seqin = self.seqinObj(summarise=False)
            ctgin = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','dna=T','seqin={0}'.format(assfile),'gapstats=F'])
            ctgdict = ctgin.seqNameDic()
            fillfas = '{0}.fillcheck.fasta'.format(basefile)
            FILLFAS = open(fillfas,'w')
            fx = 0  # Gap fills
            nx = 0  # Unchanged sequences
            ex = 0  # Edited sequences
            for seq in seqin.seqs():
                self.progLog('\r#FILL','Gap-filling sequences -> {2}: {0} unchanged; {1} edited.'.format(nx,ex,fillfas))
                (seqname, sequence) = seqin.getSeq(seq)
                sname = seqin.shortName(seq)
                if sname not in gdb.index('seqname'):
                    nx += 1
                    FILLFAS.write('>{0}\n{1}\n'.format(seqname,sequence))
                    continue
                seqlen = len(sequence)
                changes = gdb.indexEntries('seqname',sname)
                changes.reverse()
                for centry in changes:
                    ctgseq = ctgin.seqSequence(ctgdict[centry['gapctg']])
                    fillseq = ctgseq[centry['ctgstart']-1:centry['ctgend']]
                    if centry['strand'] == '-': fillseq = rje_sequence.reverseComplement(fillseq)
                    sequence = sequence[:centry['start']] + fillseq + sequence[centry['end']:]
                    fx += 1
                ex += 1
                newname = '{0}-{1}fill'.format(sname,len(changes))
                for centry in changes: centry['newname'] = newname
                self.printLog('\r#FILL','Gap-filled {0} -> {1}: {2} -> {3}.  '.format(sname,newname,rje_seqlist.dnaLen(seqlen),rje_seqlist.dnaLen(len(sequence))))
                FILLFAS.write('>{0} {1} (Diploidocus gap-filled)\n{2}\n'.format(newname,seqin.seqDesc(seq),sequence))
            FILLFAS.close()
            self.printLog('\r#FILL','Gap-filled sequences output to {2}: {0} unchanged; {1} edited; {3} of {4} gaps filled.'.format(nx,ex,fillfas,fx,cdb.entryNum()))
            gdb.saveToFile()
            if fx:
                rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','seqin={0}'.format(fillfas),'summarise=T','dna=T'])

            ## ~ [4f] Check the newly-filled gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # regcheck=FILE   : File of SeqName, Start, End positions for read coverage checking [None]
            # checkfields=LIST: Fields in checkpos file to give Locus, Start and End for checking [Hit,SbjStart,SbjEnd]
            # checkflanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
            # spanid=X        : Generate sets of read IDs that span veccheck/regcheck regions, grouped by values of field X []
            # regcnv=T/F      : Whether to calculate mean depth and predicted CNV of regcheck regions based on SCdepth [True]
            # gfftype=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
            # subforks=INT    : Number of forks for assembly subproccesses during gapfill and gapass modes [1]
            #!# Replace this with a new *.fillcheck.tdt file that has full-length filled gaps and just the edges.
            cdb = db.copyTable(gdb,'fillcheck',replace=True,add=True)
            cdb.newKey(['newname','newstart','newend'])
            for centry in cdb.entries():
                newentry = {'newstart':centry['newstart']-1,'newend':centry['newstart'],'type':'start'}
                newentry = rje.combineDict(newentry,centry,overwrite=False)
                cdb.addEntry(newentry)
                newentry = {'newstart':centry['newend'],'newend':centry['newend']+1,'type':'end'}
                newentry = rje.combineDict(newentry,centry,overwrite=False)
                cdb.addEntry(newentry)
            cdb.saveToFile()
            checkfile = '{0}.fillcheck.tdt'.format(basefile)
            checkcmd = self.cmd_list+['seqin={0}'.format(fillfas),'regcheck={0}'.format(checkfile),'checkfields=newname,newstart,newend','runmode=regcheck','paf=None','basefile={0}.fillcheck'.format(basefile)]
            if not self.getNum('SCDepth'):
                checkcmd.append('regcnv=F')
                self.printLog('#CHECK','Single-copy read depth (scdepth=NUM) not given: setting regcnv=F')
            if self.getInt('SubForks') > self.getInt('Forks'): checkcmd.append('forks={0}'.format(self.getInt('SubForks')))
            spanchecker = Diploidocus(self.log,checkcmd)
            spanchecker.run()

            ## ~ [4g] Fix any gaps without start and end support ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add reading in of *.fillcheck.checkpos.tdt and reversing gap-filling of gaps without edge support
            checkfile = '{0}.fillcheck.checkpos.tdt'.format(basefile)
            fdb = db.addTable(checkfile,mainkeys=['newname','newstart','newend','type'],name='fillchecked',ignore=[],expect=True)
            fformats = {}
            for field in fdb.fields():
                if field not in ['seqname','gapctg','strand','type','newname']: fformats[field] = 'int'
            fdb.dataFormat(fformats)
            fdb.dropEntriesDirect('type',['start','end'],inverse=True)
            fdb.dropEntriesDirect('Span0',[0],inverse=True)
            fdb.compress(['newname','newstart','newend'],default='max',rules={'type':'list'},joinchar='&')
            self.printLog('#REVGAP','{0} filled gaps identified for reversion based on inserted chunk ends without spanning reads'.format(fdb.entryNum()))
            if fdb.entryNum():     #?# Add toggle to switch this off?
                #i# Load sequences for reversion
                origin = seqin
                origdict = origin.seqNameDic()
                seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','dna=T','seqin={0}'.format(fillfas),'gapstats=F'])
                #i# Revert
                fillfas = '{0}.gapfill.fasta'.format(basefile)
                FILLFAS = open(fillfas,'w')
                rx = 0  # Reverted gap fills
                nx = 0  # Unchanged sequences
                ex = 0  # Edited sequences
                for seq in seqin.seqs():
                    self.progLog('\r#FILL','Reversing unsupported gap-filling -> {2}: {0} unchanged; {1} edited.'.format(nx,ex,fillfas))
                    (seqname, sequence) = seqin.getSeq(seq)
                    sname = seqin.shortName(seq)
                    if sname not in fdb.index('newname'):
                        nx += 1
                        FILLFAS.write('>{0}\n{1}\n'.format(seqname,sequence))
                        continue
                    seqlen = len(sequence)
                    changes = fdb.indexEntries('newname',sname)
                    changes.reverse()
                    for centry in changes:
                        ctgseq = origin.seqSequence(origdict[centry['seqname']])
                        fillseq = ctgseq[centry['start']-1:centry['end']]
                        sequence = sequence[:centry['newstart']] + fillseq + sequence[centry['newend']:]
                        rx += 1
                    ex += 1
                    newname = '{0}-{1}rev'.format(sname,len(changes))
                    #for centry in changes: centry['newname'] = newname
                    if rx: self.printLog('\r#FILL','Corrected gap-filled {0} -> {1}: {2} -> {3}.  '.format(sname,newname,rje_seqlist.dnaLen(seqlen),rje_seqlist.dnaLen(len(sequence))))
                    FILLFAS.write('>{0} {1} (Diploidocus gap-filled)\n{2}\n'.format(newname,seqin.seqDesc(seq),sequence))
                FILLFAS.close()
                self.printLog('\r#FASOUT','Corrected Gap-filled sequences output to {2}: {0} unchanged; {1} edited; {3} of {4} filled gaps reinstated.'.format(nx,ex,fillfas,rx,fx))
            else:
                checkfas = fillfas
                fillfas = '{0}.gapfill.fasta'.format(basefile)
                os.rename(checkfas,fillfas)
                self.printLog('#FASOUT','{0} renamed {1}'.format(checkfas,fillfas))
            rje_seqlist.SeqList(self.log, self.cmd_list + ['seqmode=file', 'autoload=T', 'seqin={0}'.format(fillfas), 'summarise=T', 'dna=T'])

        return True
    except:
        self.errorLog(self.zen())
    return False
#########################################################################################################################
def regionSpan(self,regfile): ### Modified version of gapspan analysis, that deals with a single region.
    '''
    ### ~ Assembly gap read-spanning analysis [runmode=gapspan] ~ ###

    This mode first identifies all the gaps in an assembly (`seqin=FILE`) (using SeqList `gapstats` or `$SEQBASE.gaps.tdt` if pre-
    existing) and then runs the read spanning analysis (`runmode=regcheck`) with `regcnv=F`. Long read data, given
    with the `reads=FILELIST` and `readtype=LIST` options, are mapped onto the assembly using minimap2 to generate a PAF file.
    This is then parsed and reads spanning each gap are identified based on their positions and the target start and end positions in the PAF file.
    In addition to absolute spanning of regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
    sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
    using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
    due to sequence length constraints.

    Spanning `spanid` output is also generated for each gap and saved in `$BASEFILE_spanid`. Each gap will be named:
    `seqname.start-end`.

    >> regfile:str = Delimited region file with seqname, start, end and (unique) name fields (i.e. like a gaps.tdt file)
    '''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        db = self.db()
        mainbase = self.baseFile()
        basefile = self.baseFile(strip_path=True)
        ## ~ [1a] ~ Check regfile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not rje.checkForFiles(filelist=[regfile],basename='',log=self.log):
            raise ValueError('Problem finding/generating region file: {0}'.format(regfile))
        checkfile = '%s.checkpos.tdt' % basefile
        ## ~ [1b] ~ Setup PAF CheckFlanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not self.list['CheckFlanks']:
            self.list['CheckFlanks'] = [0]
            self.cmd_list.append('checkflanks=0')
        minflank = self.list['CheckFlanks'][0]
        minspan = 'Span%d' % minflank
        maxflank = self.list['CheckFlanks'][-1]

        ### ~ [2] ~ Load regions and run checkpos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cdb = None
        if rje.checkForFiles(filelist=['.checkpos.tdt'],basename=basefile,log=self.log) and not self.force():
            cdb = db.addTable(checkfile,mainkeys=['seqname','start','end'],name='checkpos',ignore=[],expect=True)
            cdb.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            if minspan not in cdb.fields() and 'Span0' in cdb.fields():
                self.printLog('#SPAN','Cannot find "{0}" field in {1}: will use Span0'.format(minspan,checkfile))
                minspan = 'Span0'
            elif minspan not in cdb.fields():
                self.printLog('#SPAN','Cannot find "{0}" field in {1}: will regenerate'.format(minspan,checkfile))
                cdb = None
        if not cdb:
            cdb = db.addTable(regfile,mainkeys=['seqname','start','end'],name='region_pos',ignore=[],expect=True)
            cdb.dataFormat({'seqlen':'int','start':'int','end':'int'})
            cdb.makeField(formula='#seqname#.#start#-#end#',fieldname='region')
            cdb.saveToFile()
            self.setStr({'RegCheck':cdb.saveToFileName()})
            self.list['CheckFields'] = ['seqname','start','end']
            checkcmd = ['checkfields=seqname,start,end','spanid=name']
            ## ~ [2a] ~ Run Checkpos and SpanID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            vfile = self.getStr('RegCheck')
            pfile = self.getPAFFile()   #baseFile() + 'checkpos.paf'
            pafcmd = self.cmd_list + ['checkpos={}'.format(vfile),'pafin={}'.format(pfile)] + checkcmd
            paf = rje_paf.PAF(self.log, pafcmd)
            cdb = paf.checkPos(save=False)
            cdb.dataFormat({'MaxFlank3':'int','seqlen':'int','end':'int'})
            for centry in cdb.entries(): centry['MaxFlank3'] = centry['seqlen'] - centry['end']
            cdb.rename('checkpos')
            cdb.saveToFile(backup=False)
            #!# Clean up region_pos file #!#

        ### ~ [3] ~ Reassemble regions (runmode=gapass) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        asslist = []    # List of assembly commands to fork out and outputs (acmd, output, assembly)
        reassemble = self.getStrLC('RunMode') in ['gapass','gapfill']
        assembler = 'flye' # self.getStr('Assembler')
        filldir = '%s_regionfill/' % basefile
        assfile = '{0}{1}.assembledregions.fasta'.format(filldir,basefile)
        if self.getStrLC('RunMode') in ['gapfill'] and rje.exists(assfile) and not self.force():
            self.printLog('#REASS','{0} found (force=F): skipping assembly'.format(assfile))
            reassemble = False
        if reassemble:
            try:
                vcheck = rje.split(os.popen('{0} --version'.format(assembler)).read())[0]
                self.printLog('#PROG','{0} version: {1}'.format(assembler,vcheck))
            except:
                raise ValueError('Assembler check error - failed to run: {0} --version'.format(assembler))
        if reassemble:
            rtype = self.list['ReadType'][0]
            if rtype in ['pacbio','pac']: rtype = 'pb'
            if rtype in ['hifi','ccs']: rtype = 'hifi'
            if rtype not in ['ont','pb','hifi']:
                self.warnLog('Read Type "%s" not recognised (pb/ont/hifi): check readtype=LIST. Will use "ont".' % rtype)
                rtype = 'ont'
            if assembler == 'flye':
                if rtype == 'ont': rtype = 'nano-raw'
                elif rtype == 'hifi': rtype = 'pacbio-hifi'
                else: rtype = 'pacbio-raw'
            bamfile = self.getBamFile()
            iddir = '%s_spanid/' % basefile
            assdir = '%s_regionassemble/' % basefile
            rje.mkDir(self,assdir)
            nonx = 0; assx = 0; failx = 0; skipx = 0
            for centry in cdb.entries():
                spanner = centry['name']
                target = '{0}{1}.{2}.assembly.fasta'.format(assdir,basefile,spanner)
                if rje.exists(target) and not self.force():
                    self.printLog('#REASS','{0} found (force=F): skipped.'.format(target))
                    skipx = 0
                    continue
                spout = '%s%s.%s.span.id' % (iddir,basefile,spanner)
                ## ~ [3a] ~ Check for read ID files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not rje.exists(spout):
                    if centry[minspan] > 0:
                        self.warnLog('No {0} file but {1}={2}. May need to delete {3} and regenerate.'.format(spout,minspan,centry[minspan],checkfile))
                    elif self.v() > 1: self.printLog('No %s file' % spout)
                    nonx += 1
                    continue
                ## ~ [3b] ~ Extract the reads from the BAM file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #!# Convert this to a forked process too?
                alog = '{0}{1}.{2}.assemble.log'.format(assdir,basefile,spanner)
                self.printLog('#BAM','Extracting region-spanning reads from BAM file.')
                # Exclude secondary and umapped reads
                region = '{0}:{1}-{2}'.format(centry['seqname'],max(centry['start']-maxflank,1),centry['end']+maxflank)
                idbam = '%s%s.%s.bam' % (assdir,basefile,spanner)
                fasout = '%s%s.%s.reads.fasta' % (assdir,basefile,spanner)
                if self.force() or not rje.exists(fasout):
                    region2bam = 'samtools view -h -F 0x104 {0} {1} | grep -e @ -f {2} | samtools view -h -b -o {3} -'.format(bamfile,region,spout,idbam)
                    logline = self.loggedSysCall(region2bam,alog,append=True,stderr=False,nologline='No stdout from samtools view',threaded=False)
                    bam2fas = 'samtools fasta {0} > {1}'.format(idbam,fasout)
                    logline = self.loggedSysCall(bam2fas,alog,append=True,stderr=False,nologline='No stdout from samtools fasta',threaded=False)
                gensize = maxflank*2
                if rje.exists(fasout) and open(fasout,'r').read():
                    fascmd = ['dna=T','raw=T','gapstat=F','summarise=F','autoload=T','seqmode=list','seqin={0}'.format(fasout)]
                    fasta = rje_seqlist.SeqList(self.log,self.cmd_list+fascmd)
                    if fasta.seqNum() < 1:
                        self.warnLog('{0} has no reads despite {1} spanning read ID output'.format(fasout,spout))
                        continue
                    fasdat = fasta.summarise(sumdb=False,save=False)
                    gensize = int(fasdat['MeanLength'])
                else:
                    self.warnLog('{0} has no reads despite {1} spanning read ID output'.format(fasout,spout))
                    continue
                ## ~ [3c] Assemble reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Adding forking
                #i# asslist = []    # List of assembly commands to fork out and outputs (acmd, output, assembly)
                if assembler == 'flye':
                    acmd = 'flye --{0} {1} --out-dir {2}{3}.{4}_flye --genome-size {5} --threads {6}'.format(rtype,fasout,assdir,basefile,spanner,gensize,self.getInt('SubForks'))
                    assembly = '{0}{1}.{2}_flye/assembly.fasta'.format(assdir,basefile,spanner)
                if self.v() > 0:
                    acmd = '{0} 2>&1 | tee -a {1}'.format(acmd,alog)
                else:
                    acmd = '{0} 2>&1 >> {1}'.format(acmd,alog)
                asslist.append((acmd,assembly,'{0}{1}.{2}.assembly.fasta'.format(assdir,basefile,spanner)))
                continue
            self.printLog('#REASS','{0} regions without read IDs'.format(rje.iStr(nonx)))
            if skipx:
                self.printLog('#REASS','{0} existing assemblies skipped (force=F).'.format(rje.iStr(skipx)))

            ### ~ [4] ~ Fork out assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#FORK','{0} region assemblies to fork out'.format(rje.iLen(asslist)))
            ## ~ [4a] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            forker = self.obj['Forker']
            forker.list['ToFork'] = []
            for afork in asslist:
                forker.list['ToFork'].append(afork[0])
            self.debug('\n'.join(forker.list['ToFork']))
            ## ~ [4b] Fork out assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if forker.list['ToFork']:
                if self.getNum('Forks') < 1:
                    #i# Warn lack of forking
                    self.printLog('#FORK','Note: program can be accelerated using forks=INT.')
                    for acmd in forker.list['ToFork']:
                        alog = afork[2][:-6] + '.log'
                        logline = self.loggedSysCall(acmd,alog,append=True,nologline='No stdout from flye',threaded=False)
                        #self.printLog('#SYS',forkcmd)
                        #os.system(forkcmd)
                else:
                    self.printLog('#FORK','Forking assemblies using {0} x {1} threads'.format(self.getNum('Forks'),self.getNum('SubForks')))
                    if forker.run():
                        self.printLog('#FORK','Forking of assemblies completed.')
                    else:
                        try:
                            self.errorLog('Assembly forking did not complete',printerror=False,quitchoice=True)
                        except:
                            raise RuntimeError('Assembly forking did not complete')
            ## ~ [4c] Process assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            assx = 0; failx = 0
            for afork in asslist:
                assembly = afork[1]
                target = afork[2]
                if rje.exists(assembly):
                    self.printLog('#REASS','{0} generated -> {1}'.format(assembly,target))
                    rje.fileTransfer(assembly,target,deletefrom=False,append=False)
                    assx += 1
                else:
                    self.printLog('#FAIL','No {0} generated.'.format(assembly))
                    failx += 1
            self.printLog('#REASS','{0} region assemblies generated; {1} failed'.format(rje.iStr(assx),rje.iStr(failx)))

        ### ~ [5] ~ Fill regions (runmode=gapfill) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #i# The idea here is to use GABLAM-style PAF mapping to map the assembled regions back onto the assembly
        if self.getStrLC('RunMode') in ['gapfill']:
            ## ~ [4a] ~ Compile the assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            assdir = '%s_regionassemble/' % basefile
            filldir = '%s_regionfill/' % basefile
            rje.mkDir(self,filldir)
            assfile = '{0}{1}.assembledregions.fasta'.format(filldir,basefile)
            if self.force() or not rje.exists(assfile):
                rje.backup(self,assfile)
                ASSFILE = open(assfile,'w')
                nonx = 0; assx = 0; nullx = 0; ctgx = 0
                for centry in cdb.entries():
                    spanner = centry['name']
                    assembly = '{0}{1}.{2}.assembly.fasta'.format(assdir,basefile,spanner)
                    if not rje.exists(assembly): nonx += 1; continue
                    aseq = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=list','autoload=F','seqin={0}'.format(assembly)])
                    aseq.loadSeq()
                    if aseq.list['Seq']:
                        sx = 1; assx += 1
                        for (name, sequence) in aseq.list['Seq']:
                            ASSFILE.write('>{0}.ctg{1} {2}\n{3}\n'.format(spanner,sx,name,sequence))
                            sx += 1; ctgx += 1
                    else: nullx += 1
                ASSFILE.close()
                self.printLog('#REGCTG','Collated {0} assembled contigs for {1} regions; {2} assemblies w/o contigs; {3} w/o assemblies'.format(rje.iStr(ctgx),rje.iStr(assx),rje.iStr(nullx),rje.iStr(nonx)))
            else:
                self.printLog('#REGCTG','Collated assembled contigs found (force=F): {0}'.format(assfile))
            ## ~ [4b] ~ Map region assemblies to genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fillbase = '%s%s.regionctg' % (filldir,self.baseFile())
            pafin = '%s.paf' % (fillbase)
            paf = rje_paf.PAF(self.log, self.cmd_list+['pafin=%s' % pafin,'%s%s.regionctg' % (filldir,self.baseFile()),'seqin=%s' % assfile,'reference=%s' % self.getStr('SeqIn')])
            defaults = paf_defaults
            paf.dict['MapOpt'] = rje.combineDict(defaults,paf.dict['MapOpt'],overwrite=True)
            paf.setup()
            if not rje.exists(pafin) or self.force():
                rje.backup(self,pafin)
                paf.minimap2()
                #!# Add localaln=T to PAF #!#
            #?# Should paf.db() basename be changed to pafin basefile for the initial output?
            ## ~ [4c] ~ GapFiller local hits table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gapfiller = '%s.regionfiller.tdt' % (fillbase)
            paf.basefile(fillbase)
            if rje.exists(gapfiller) and not self.force():
                udb = paf.db().addTable(gapfiller,mainkeys=['Hit','SbjStart','SbjEnd'],name='regionfiller',ignore=[],expect=True)
                udb.dataFormat({'QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','QryLen':'int','SbjLen':'int'})
            else:
                paf.run()
                paf.db('hitunique').rename('regionfiller')     # ['Qry','QryStart','QryEnd']
                paf.db('regionfiller').dropFields(['nn','cs'])
                for entry in paf.db('regionfiller').entries():
                    if entry['SbjStart'] > entry['SbjEnd']:
                        (entry['SbjStart'],entry['SbjEnd']) = (entry['SbjEnd'],entry['SbjStart'])
                paf.db('regionfiller').newKey(['Hit','SbjStart','SbjEnd'])
                paf.db('regionfiller').saveToFile()
                #  Qry     Hit     AlnNum  BitScore        Expect  Length  Identity        Positives       QryStart        QryEnd  SbjStart        SbjEnd  QryLen  SbjLen  Strand
            self.printLog('#REGMAP','%s%s.regionctg.* output generated.' % (filldir,self.baseFile()))
            ## ~ [4d] ~ Identify region fills using local hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Work through by qry position, looking for a single hit or pair of hits that spans the region
            #i# -> For a single hit, replace the whole local hit region
            #i# -> For a pair of hits, replace the region between the Sbj regions
            #i# => Build a list of (seqname, start, end, regionctg, start, end, strand)
            #i# .. Reverse sort of replace sequences in turn
            #i# => Output to *.regionfill.fasta
            #i# -> Report number of regions filled and new stats
            paf.basefile(mainbase)
            udb = paf.db('regionfiller')
            preventry = None
            regionfills = []   # Tuples of (seqname, start, end, regionctg, start, end, strand)
            filled = []
            #i# Filter on local alignment length (and identity?)
            minloclen = self.getInt('MinLocLen')
            minlocid = self.getPerc('MinLocID')
            if minloclen > 1: udb.dropEntries('Length<%d' % minloclen)
            if minlocid > 0:
                udb.makeField('100.0*Identity/Length','LocID')
                udb.dropEntries('LocID<%s' % minlocid)
                udb.dropField('LocID')
            #i# Identify regions to fill
            ex = 0.0; etot = udb.entryNum()
            for entry in udb.entries(sorted=True):
                self.progLog('\r#FILL','Identifying region-filling re-assemblies: %.2f%%' % (ex/etot)); ex += 100.0
                regionid = '.'.join(entry['Qry'].split('.')[:-1])  # PTEXCHR1.01.10697073-10697172.ctg1
                if regionid in filled: continue
                #i# Check for matching sequences
                if not entry['Qry'].startswith('{0}.'.format(entry['Hit'])): continue
                #i# Check preventry
                if preventry and preventry['Qry'] != entry['Qry']: preventry = None
                #i# Check for single hit spanning region
                #if entry['Strand'] == '-': (entry['SbjStart'],entry['SbjEnd']) = (entry['SbjEnd'],entry['SbjStart'])
                [regionx,regiony] = entry['Qry'].split('.')[-2].split('-')
                regionx = int(regionx)
                regiony = int(regiony)
                if entry['SbjStart'] < regionx and entry['SbjEnd'] > regiony:
                    filled.append(regionid)
                    regionfills.append((entry['Hit'],entry['SbjStart'],entry['SbjEnd'],entry['Qry'],entry['QryStart'],entry['QryEnd'],entry['Strand'],'span'))
                    preventry = entry
                    continue
                if not preventry:
                    preventry = entry
                    continue
                #i# Check for dual hit spanning region
                if preventry['SbjEnd'] <= regionx and entry['SbjStart'] >= regiony and preventry['Strand'] == entry['Strand']:
                    if (entry['Strand'] == '+' and entry['SbjStart'] > preventry['SbjEnd']):
                        filled.append(regionid)
                        regionfills.append((entry['Hit'],preventry['SbjEnd'],entry['SbjStart'],entry['Qry'],preventry['QryEnd'],entry['QryStart'],entry['Strand'],'flank'))
                    elif (entry['Strand'] == '-' and preventry['SbjStart'] > entry['SbjEnd']):
                        filled.append(regionid)
                        regionfills.append((entry['Hit'],preventry['SbjEnd'],entry['SbjStart'],entry['Qry'],entry['QryEnd'],preventry['QryStart'],entry['Strand'],'flank'))
                preventry = entry
            self.progLog('\r#FILL','Identified %s region-filling re-assemblies from %s.' % (rje.iLen(regionfills),regionfiller))
            ## ~ [4e] ~ RegionFill using local hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add a check for overlapping replacements? #!#
            regionfills.sort()
            regionshift = {}
            gdb = db.addEmptyTable('regionfill',['seqname','start','end','regionctg','ctgstart','ctgend','strand','regionlen','type','newname','newstart','newend','newchunklen'],['seqname','start','end'],log=False)
            for filldat in regionfills:
                fentry = {'seqname':filldat[0],'start':filldat[1],'end':filldat[2],'regionctg':filldat[3],'ctgstart':filldat[4],'ctgend':filldat[5],'strand':filldat[6],'type':filldat[7]}
                if fentry['seqname'] not in regionshift: regionshift[fentry['seqname']] = 0
                fentry['newname'] = fentry['seqname']
                fentry['newstart'] = fentry['start'] + regionshift[fentry['seqname']]
                regionshift[fentry['seqname']] += (fentry['ctgend'] - fentry['ctgstart']) - (fentry['end'] - fentry['start'])
                fentry['newend'] = fentry['end'] + regionshift[fentry['seqname']]
                fentry['newchunklen'] = fentry['newend'] - fentry['newstart'] + 1
                [regionx,regiony] = fentry['regionctg'].split('.')[-2].split('-')
                fentry['regionlen'] = int(regiony) - int(regionx) + 1
                gdb.dict['Data'][(filldat[0],filldat[1],filldat[2])] = fentry
            gdb.index('seqname')
            seqin = self.seqinObj(summarise=False)
            ctgin = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','dna=T','seqin={0}'.format(assfile),'gapstats=F'])
            ctgdict = ctgin.seqNameDic()
            fillfas = '{0}.fillcheck.fasta'.format(basefile)
            FILLFAS = open(fillfas,'w')
            fx = 0  # region fills
            nx = 0  # Unchanged sequences
            ex = 0  # Edited sequences
            for seq in seqin.seqs():
                self.progLog('\r#FILL','Region-filling sequences -> {2}: {0} unchanged; {1} edited.'.format(nx,ex,fillfas))
                (seqname, sequence) = seqin.getSeq(seq)
                sname = seqin.shortName(seq)
                if sname not in gdb.index('seqname'):
                    nx += 1
                    FILLFAS.write('>{0}\n{1}\n'.format(seqname,sequence))
                    continue
                seqlen = len(sequence)
                changes = gdb.indexEntries('seqname',sname)
                changes.reverse()
                for centry in changes:
                    ctgseq = ctgin.seqSequence(ctgdict[centry['regionctg']])
                    fillseq = ctgseq[centry['ctgstart']-1:centry['ctgend']]
                    if centry['strand'] == '-': fillseq = rje_sequence.reverseComplement(fillseq)
                    sequence = sequence[:centry['start']] + fillseq + sequence[centry['end']:]
                    fx += 1
                ex += 1
                newname = '{0}-{1}fill'.format(sname,len(changes))
                for centry in changes: centry['newname'] = newname
                self.printLog('\r#FILL','Region-filled {0} -> {1}: {2} -> {3}.  '.format(sname,newname,rje_seqlist.dnaLen(seqlen),rje_seqlist.dnaLen(len(sequence))))
                FILLFAS.write('>{0} {1} (Diploidocus region-filled)\n{2}\n'.format(newname,seqin.seqDesc(seq),sequence))
            FILLFAS.close()
            self.printLog('\r#FILL','Region-filled sequences output to {2}: {0} unchanged; {1} edited; {3} of {4} regions filled.'.format(nx,ex,fillfas,fx,cdb.entryNum()))
            gdb.saveToFile()
            if fx:
                rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','seqin={0}'.format(fillfas),'summarise=T','dna=T'])

            ## ~ [4f] Check the newly-filled regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # regcheck=FILE   : File of SeqName, Start, End positions for read coverage checking [None]
            # checkfields=LIST: Fields in checkpos file to give Locus, Start and End for checking [Hit,SbjStart,SbjEnd]
            # checkflanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
            # spanid=X        : Generate sets of read IDs that span veccheck/regcheck regions, grouped by values of field X []
            # regcnv=T/F      : Whether to calculate mean depth and predicted CNV of regcheck regions based on SCdepth [True]
            # gfftype=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
            # subforks=INT    : Number of forks for assembly subproccesses during regionfill and regionass modes [1]
            #!# Replace this with a new *.fillcheck.tdt file that has full-length filled regions and just the edges.
            cdb = db.copyTable(gdb,'fillcheck',replace=True,add=True)
            cdb.newKey(['newname','newstart','newend'])
            for centry in cdb.entries():
                newentry = {'newstart':centry['newstart']-1,'newend':centry['newstart'],'type':'start'}
                newentry = rje.combineDict(newentry,centry,overwrite=False)
                cdb.addEntry(newentry)
                newentry = {'newstart':centry['newend'],'newend':centry['newend']+1,'type':'end'}
                newentry = rje.combineDict(newentry,centry,overwrite=False)
                cdb.addEntry(newentry)
            cdb.saveToFile()
            checkfile = '{0}.fillcheck.tdt'.format(basefile)
            checkcmd = self.cmd_list+['seqin={0}'.format(fillfas),'regcheck={0}'.format(checkfile),'checkfields=newname,newstart,newend','runmode=regcheck','paf=None','basefile={0}.fillcheck'.format(basefile)]
            if not self.getNum('SCDepth'):
                checkcmd.append('regcnv=F')
                self.printLog('#CHECK','Single-copy read depth (scdepth=NUM) not given: setting regcnv=F')
            if self.getInt('SubForks') > self.getInt('Forks'): checkcmd.append('forks={0}'.format(self.getInt('SubForks')))
            spanchecker = Diploidocus(self.log,checkcmd)
            spanchecker.run()

            ## ~ [4g] Fix any regions without start and end support ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add reading in of *.fillcheck.checkpos.tdt and reversing region-filling of regions without edge support
            checkfile = '{0}.fillcheck.checkpos.tdt'.format(basefile)
            fdb = db.addTable(checkfile,mainkeys=['newname','newstart','newend','type'],name='fillchecked',ignore=[],expect=True)
            fformats = {}
            for field in fdb.fields():
                if field not in ['seqname','regionctg','strand','type','newname']: fformats[field] = 'int'
            fdb.dataFormat(fformats)
            fdb.dropEntriesDirect('type',['start','end'],inverse=True)
            fdb.dropEntriesDirect('Span0',[0],inverse=True)
            fdb.compress(['newname','newstart','newend'],default='max',rules={'type':'list'},joinchar='&')
            self.printLog('#REVREG','{0} filled regions identified for reversion based on inserted chunk ends without spanning reads'.format(fdb.entryNum()))
            if fdb.entryNum():     #?# Add toggle to switch this off?
                #i# Load sequences for reversion
                origin = seqin
                origdict = origin.seqNameDic()
                seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','dna=T','seqin={0}'.format(fillfas),'gapstats=F'])
                #i# Revert
                fillfas = '{0}.regionfill.fasta'.format(basefile)
                FILLFAS = open(fillfas,'w')
                rx = 0  # Reverted region fills
                nx = 0  # Unchanged sequences
                ex = 0  # Edited sequences
                for seq in seqin.seqs():
                    self.progLog('\r#FILL','Reversing unsupported region-filling -> {2}: {0} unchanged; {1} edited.'.format(nx,ex,fillfas))
                    (seqname, sequence) = seqin.getSeq(seq)
                    sname = seqin.shortName(seq)
                    if sname not in fdb.index('newname'):
                        nx += 1
                        FILLFAS.write('>{0}\n{1}\n'.format(seqname,sequence))
                        continue
                    seqlen = len(sequence)
                    changes = fdb.indexEntries('newname',sname)
                    changes.reverse()
                    for centry in changes:
                        ctgseq = origin.seqSequence(origdict[centry['seqname']])
                        fillseq = ctgseq[centry['start']-1:centry['end']]
                        sequence = sequence[:centry['newstart']] + fillseq + sequence[centry['newend']:]
                        rx += 1
                    ex += 1
                    newname = '{0}-{1}rev'.format(sname,len(changes))
                    #for centry in changes: centry['newname'] = newname
                    if rx: self.printLog('\r#FILL','Corrected region-filled {0} -> {1}: {2} -> {3}.  '.format(sname,newname,rje_seqlist.dnaLen(seqlen),rje_seqlist.dnaLen(len(sequence))))
                    FILLFAS.write('>{0} {1} (Diploidocus region-filled)\n{2}\n'.format(newname,seqin.seqDesc(seq),sequence))
                FILLFAS.close()
                self.printLog('\r#FASOUT','Corrected region-filled sequences output to {2}: {0} unchanged; {1} edited; {3} of {4} filled regions reinstated.'.format(nx,ex,fillfas,rx,fx))
            else:
                checkfas = fillfas
                fillfas = '{0}.regionfill.fasta'.format(basefile)
                os.rename(checkfas,fillfas)
                self.printLog('#FASOUT','{0} renamed {1}'.format(checkfas,fillfas))
            rje_seqlist.SeqList(self.log, self.cmd_list + ['seqmode=file', 'autoload=T', 'seqin={0}'.format(fillfas), 'summarise=T', 'dna=T'])

        return True
    except:
        self.errorLog(self.zen())
    return False
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: 'MAIN' PROGRAM                                                                                         #
#########################################################################################################################
if __name__ == "__main__":  ### Print message to screen if called from commandline.
    try: print( __doc__)
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
