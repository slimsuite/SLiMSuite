#!/usr/bin/python

# See below for name and description
# Copyright (C) 2010 Richard J. Edwards <software@cabbagesofdoom.co.uk>
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
# Author contact: <software@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       RJE_HTML
Description:  Module for generating HTML 
Version:      0.3.0
Last Edit:    23/04/18
Copyright (C) 2010  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is primarily for general methods for making HTML pages for other modules. 

Commandline:
    stylesheets=LIST    : List of CSS files to use ['http://www.slimsuite.unsw.edu.au/stylesheets/slimhtml.css']
    tabber=FILE         : Tabber javascript file location ['tabber.js']
    border=X            : Border setting for tables [0]
    nobots=T/F          : Whether to add code to avoid Google Bots [True]
    analytics=X         : Google Analytics code to use with pages []
    javascript=PATH     : Path to javascript files for tabs etc. ['http://www.slimsuite.unsw.edu.au/javascript/']
    jscripts=LIST       : List of javascript files to load ['stupidtable.js?dev']
    keywords=LIST       : List of keywords to put in page metadata []
    title=X             : Default title for HTML page []
    copyright=X         : Copyright statement for page ['RJ Edwards 2015']

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_slim, rje_uniprot, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added additional commandline options, including Google Analytics
    # 0.2.0 - Added delimited text to HTML table conversion.
    # 0.2.1 - Updated default CSS to http://www.slimsuite.unsw.edu.au/stylesheets/slimhtml.css.
    # 0.3.0 - Added optional loading of javascript files and stupidtable.js?dev default.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add Google Analytics
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, cyear) = ('RJE_HTML', '0.3.0', 'April 2018', '2010')
    description = 'RJE HTML Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program when called from commandline.
    '''
    Basic Setup of Program when called from commandline:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:### ~ [1] ~ Initial Command Setup & Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = makeInfo()                                   # Sets up Info object with program details
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2 
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Problem during initial setup.'; raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: HTML Class                                                                                              #
#########################################################################################################################
class HTML(rje.RJE_Object):     
    '''
    HTML Class. Author: Rich Edwards (2010).

    Info:str
    - Analytics = Google Analytics code to use with pages []
    - Copyright = Copyright statement for page ['RJ Edwards 2012']
    - Javascript = Path to javascript files for tabs etc. ['http://www.slimsuite.unsw.edu.au/javascript/']
    - Tabber = Tabber javascript file location ['tabber.js']
    - Title = Default title for HTML page []
    
    Opt:boolean
    - NoBots = Whether to avoid Google Bot trawlers [True]

    Stat:numeric
    - Border = Border setting for tables [0]

    List:list
    - JScripts=LIST       : List of javascript files to load ['stupidtable.js?dev']
    - Keywords = List of keywords to put in page metadata []
    - StyleSheets = List of CSS files to use [http://www.slimsuite.unsw.edu.au/stylesheets/slimhtml.css]

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Analytics','Copyright','Javascript','Tabber','Title']
        self.optlist = ['NoBots']
        self.statlist = ['Border']
        self.listlist = ['JScripts','Keywords','StyleSheets']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.HTMLdefaults()
    def HTMLdefaults(self):
        self.setInfo({'Javascript':'http://www.slimsuite.unsw.edu.au/javascript/','Tabber':'../tabber.js',
                      'Copyright':'RJ Edwards 2015'})
        self.setStat({'Border':0})
        self.setOpt({'NoBots':True})
        self.list['StyleSheets'] = ['http://www.slimsuite.unsw.edu.au/stylesheets/slimhtml.css']
        self.list['JScripts'] = ['stupidtable.js?dev']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#########################################################################################################################
    def _cmdList(self):  self.HTMLcmdList()    ### Sets Attributes from commandline
    def HTMLcmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdReadList(cmd,'info',['Analytics','Copyright','Javascript','Tabber','Title'])
                self._cmdReadList(cmd,'int',['Border'])
                self._cmdReadList(cmd,'list',['JScripts','StyleSheets'])
                self._cmdReadList(cmd,'opt',['NoBots'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.info['Analytics'].lower() in ['','none']: self.info['Analytics'] = ''
        if self.info['Copyright'].lower() in ['','none']: self.info['Copyright'] = ''
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def tabberHTML(self,id,tablist,level=0):     ### Returns text for Tabber HTML
        '''Returns text for Tabber HTML.'''
        try: return tabberHTML(id,tablist,level)
        except: self.errorLog('TabberHTML Error')
        return '!Error!'
#########################################################################################################################
    def htmlHead(self,title=None,tabber=True,frontpage=False,keywords=[],redirect='',refresh=0):    ### Returns text for top of HTML file
        #!# Add additional commandline options for keywords, title, javascript etc.
        if not title: title = self.getStr('Title')
        if not keywords: keywords = self.list['Keywords']
        return htmlHead(title,self.list['StyleSheets'],tabber,frontpage,self.getOpt('NoBots'),keywords,self.info['Javascript'],self.info['Analytics'],redirect,refresh,self.list['JScripts'])
#########################################################################################################################
    def htmlTail(self,tabber=True,stupidtable=False):  ### Returns text for bottom of HTML
        return htmlTail(self.getStr('Copyright'),tabber,stupidtable)
#########################################################################################################################
### End of SECTION II: HTML Class                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def htmlHead(title,stylesheets=['../example.css','../redwards.css'],tabber=True,frontpage=False,nobots=True,keywords=[],javascript='../javascript/',analytics='',redirect='',refresh=0,jscripts=['stupidtable.js?dev']):    ### Returns text for top of HTML file
    '''
    Returns text for top of HTML file.
    >> title:str = Title of webpage.
    >> stylesheets:list = List of stylesheets to use.
    >> tabber:bool [True] = whether page has tabber tabs
    >> frontpage:bool [False] = whether to replace all '../' links with './'
    >> nobots:bool [True] = whether to screen page from bot discovery
    >> keywords:list [] = List of keywords for header
    >> javascript:str ['../javascript/'] = Path to javascript files for tabber tabs
    >> redirect:str [''] = URL to redirect to.
    >> refresh:int [0] = Time to redirect
    >> jscripts:list ['stupidtable.js?dev'] = List of javascript files to load.
    '''
    html = ['<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">',
            '<html lang="en">','','<!-- ~~~~~~~~~~~~~~~ HTML head data ~~~~~~~~~~~~~~~~~ -->','<head>',
            '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">', #!# Add additional metadata #!#
            '<title>%s</title>' % title,'']
    if keywords: html.insert(6,'<meta http-equiv="Keywords" name="Keywords" content="%s">' % string.join(keywords, ', '))
    if nobots: html.insert(6,'<meta name="robots" content="index,nofollow">')
    for stylesheet in stylesheets:
        if frontpage: html.append('<link rel="stylesheet" href="%s" TYPE="text/css" MEDIA="screen">' % string.replace(stylesheet,'../','./'))
        else: html.append('<link rel="stylesheet" href="%s" TYPE="text/css" MEDIA="screen">' % stylesheet)
    if tabber:
        html += ['','<!-- ~~~~~~~~~~~ Tabber Javascript ~~~~~~~~~~~ -->','<script type="text/javascript">',
                 'document.write(\'<style type="text/css">.tabber{display:none;}<\/style>\');',
                 'var tabberOptions = {',' /* Optional: start manually (run tabber) at end of file','*/',
                 '\'manualStartup\':true','};','</script>','<!-- Load the tabber code -->']
        if frontpage and javascript == '../javascript/': html += ['<script type="text/javascript" src="./javascript/tabber.js"></script>','']
        else: html += ['<script type="text/javascript" src="%stabber.js"></script>' % javascript,'']
    if analytics:
        html += ['','<!-- ~~~~~~~~~~~ Google Analytics Script ~~~~~~~~~ -->','<script type="text/javascript">',
                 'var _gaq = _gaq || [];',"_gaq.push(['_setAccount', '%s']);" % analytics,"_gaq.push(['_trackPageview']);",
                 '','(function() {',"  var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;",
                 "  ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';",
                 "  var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);",'})();','</script>']
    if jscripts:
        html += ['','<!-- ~~~~~~~~~~~ Google Analytics Script ~~~~~~~~~ -->']
        for js in jscripts: html += ['<script src="%s%s"></script>' % (javascript,js)]
    if redirect: html.append(redirectToURL(redirect,refresh))
    html += ['</head>','<!-- ~~~~~~~~~~~~~~~ End of HTML head data ~~~~~~~~~~~~~~~~~ -->','','<body>','']
    #print string.join(html,'\n')
    return string.join(html,'\n')
#########################################################################################################################
def htmlTail(copyright='RJ Edwards 2012',tabber=True,stupidtable=False):  ### Returns text for bottom of HTML
    '''
    Returns text for bottom of HTML.
    >> copyright:str = copyright text'
    >> tabber:bool = whether page has tabber tabs
    '''
    t = string.split(time.asctime(time.localtime(time.time())))
    datetime = '%s %s %s' % (t[2],t[1],t[-1])
    html = ['','<!-- ~~~~~~~~~~~~~~ HTML tail data ~~~~~~~~~~~~~~~~~ -->',
            '<HR><FONT COLOR=#979E45 SIZE=2>&copy; %s. Last modified %s.</FONT></P>' % (copyright,datetime),'']
    if tabber:
        html += ['<script type="text/javascript">','/* manualStartup=true so need to run it now */',
                 'tabberAutomatic(tabberOptions);','</script>','']
    #!# Need to replace this with custom script options
    if stupidtable:
        html += ['<script>','$(".instances").stupidtable();','</script>','']
    html += ['</body>','</html>','<!-- ~~~~~~~~~~~~~~ End of HTML tail data ~~~~~~~~~~~~~~~~~ -->']
    #print string.join(html,'\n')
    return string.join(html,'\n')
#########################################################################################################################
def tabberHTML(tabid,tablist,level=0):     ### Returns text for Tabber HTML
    '''
    Returns text for Tabber HTML.
    >> tabid:str = Identifier for Tabber object
    >> tablist:list = List of (tab_id, tab_html_text[, tab_title]) tuples
    >> level:int = Level of Tabber object (base = level)
    '''
    jointxt = '\n' + '    ' * level 
    html = ['<!-- ~~~~~~~~~~~~~~~ %s Tabber Div ~~~~~~~~~~~~~~~ -->' % tabid,'','<div class="tabber" id="%s">' % tabid,'']
    #print html
    for tab in tablist:
        #print tab
        #print tab[0],tab[1]
        #print tabberTabHTML(tab[0],tab[1])
        if len(tab) > 2: html += string.split(tabberTabHTML(tab[0],tab[1],tab[2]),'\n')
        else: html += string.split(tabberTabHTML(tab[0],tab[1]),'\n')
    html += ['</div>','<!-- ~~~~~~~~~~~~~~~ End of %s Tabber Div ~~~~~~~~~~~~~~~ -->' % tabid,]
    ### Join HTML code ###
    htmljoin = ''
    pre = 0
    while html:
        nextline = html.pop(0)
        if pre <= 0: htmljoin += '    ' * level; pre = 0
        htmljoin += nextline + '\n'
        pre += nextline.lower().count('<pre>')
        pre -= nextline.lower().count('</pre>')
    return htmljoin
#########################################################################################################################
def tabberTabHTML(tabid,text,title=''):          ### Returns text for TabberTab HTML
    '''
    Returns text for TabberTab HTML.
    >> title:str = Text for title of TabberTab
    >> text:str = HTML text for TabberTab content
    '''
    if not title: title = tabid
    html = ['','<!-- ~~~ %s TabberTab div ~~~ -->' % tabid,'<div class="tabbertab" title="%s" id="%s">' % (title,tabid),'']
    html += string.split(text,'\n')
    html += ['','</div>','<!-- ~~~ %s TabberTab end ~~~ -->' % tabid]
    #print string.join(html,'\n  ')
    if string.join(html).upper().find('<PRE>') >= 0: return string.join(html,'\n')       
    else: return string.join(html,'\n  ')
#########################################################################################################################
def geneLink(gene,frontpage=False):     ### Returns gene link text
    '''Returns gene link text.'''
    if frontpage: return '<a href="./gene/%s.html">%s</a>' % (gene,gene)
    else: return '<a href="../gene/%s.html">%s</a>' % (gene,gene)
#########################################################################################################################
def domainLink(domain,frontpage=False):     ### Returns gene link text
    '''Returns domain link text.'''
    if frontpage: return '<a href="./domain/%s.html">%s</a>' % (domain,domain)
    else: return '<a href="../domain/%s.html">%s</a>' % (domain,domain)
#########################################################################################################################
def slimLink(pattern,frontpage=False):  ### Returns SLiM link text
    '''Returns gene link text.'''
    if frontpage: return '<a href="./slim/%s.html">%s</a>' % (rje_slim.slimFromPattern(pattern),pattern)
    else: return '<a href="../slim/%s.html">%s</a>' % (rje_slim.slimFromPattern(pattern),pattern)
#########################################################################################################################
def seqDetailsHTML(callobj,gene,dbxref):    #gene,seqid,dbxref,desc,godata):  ### Returns HTML text for seq details table
    '''
    Returns HTML text for seq details table.
    >> gene:str = Gene symbol
    >> seqid:str = Sequence Identifier
    >> dbxref:dict = Dictionary of {db:id} for GeneCards, EBI, EnsEMBL, HPRD, OMIM
    >> desc:str = Sequence description
    >> godata:dict = {CC/BP/MF:[(id,name)] list}
    '''
    try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqid = '-'; desc = '<i>No sequence mapping.</i>'; godata = {}
        if dbxref:      #gene in callobj.data('DBXRef'):
            #dbxref = callobj.data('DBXRef')[gene]
            seqid = dbxref['EnsLoci']
            desc = dbxref['EnsDesc']
            ensg = dbxref['EnsEMBL']
            try: godata = callobj.dict['GO'][ensg]
            except: pass

        gfont = '<FONT SIZE=6 FACE="Verdana" COLOR=#014359>'
        ifont = '<FONT SIZE=5 FACE="Verdana" COLOR=#979E45>'
        dfont = '<FONT SIZE=4 COLOR=#014359>'
        html = ['','<!-- ~~~~ %s Gene Summary details table ~~~~ -->' % gene,'<table width="100%">',
                '<tr valign="top">'  #<td width="80%">', '<table>','<tr>',
                '<td width="30%%"><a href="../gene/%s.html">%s<b>%s</b></FONT></a></td>' % (gene,gfont,gene),
                '<td width="50%%">%s<b>%s</b></FONT></td>' % (ifont,seqid),
                '<td width="20%" align="right" rowspan="3">',
                '<a href="../index.htm"><img src="../resources/SBS_100.png" height="100" align="RIGHT" border="0" alt="Home"></a>',
                '</td></tr>']
        #x#if 'Gene' not in dbxref: dbxref['Gene'] = gene
        html += ['<tr><td colspan=2>']
        for db in ['Gene','UniProt','EnsEMBL','HPRD','OMIM']:
            if db not in dbxref: continue
            if db == 'Gene':
                href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' % dbxref[db]
                alt = 'GeneCards %s' % dbxref[db]
                src = '../resources/genecard.png'
            if db == 'UniProt':
                href = 'http://www.uniprot.org/uniprot/%s' % dbxref[db]
                alt = 'UniProt %s' % dbxref[db]
                src = '../resources/logo_ebi.png'
            if db == 'EnsEMBL':
                href = 'http://www.ensembl.org/Homo_sapiens/geneview?gene=%s' % dbxref[db]
                alt = 'EnsEMBL %s' % dbxref[db]
                src = '../resources/e-bang.gif'
            if db == 'HPRD':
                href = 'http://www.hprd.org/summary?protein=%s&isoform_id=%s_1&isoform_name=Isoform_1' % (dbxref[db],dbxref[db])
                alt = 'HPRD %s' % dbxref[db]
                src = '../resources/hprd.png'
            if db == 'OMIM':
                href = 'http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s' % dbxref[db]
                alt = 'OMIM %s' % dbxref[db]
                src = '../resources/omim.png'
            html += ['<a href="%s"><img alt="%s" src="%s" align="BOTTOM" border="0" height="50"></a>' % (href,alt,src)] 
        html += ['</td></tr>','<tr><td colspan=2>%s<p>%s</p></FONT></td></tr>' % (dfont,desc)]
        html += ['</table>','<!-- ~~~~ End %s Gene Summary details table ~~~~ -->' % gene,'']
        gtab = []
        for gtype in ['CC','BP','MF']:
            if gtype in godata:
                gdict = {}
                for gotup in godata[gtype]:
                    if gotup[1] in ['cellular_component','biological_process','molecular_function']: continue
                    gdict[gotup[1]] = '<a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s">%s</a>' % (gotup[0],gotup[1])
                if gdict:
                    ghtml = []
                    for g in rje.sortKeys(gdict): ghtml.append(gdict[g])
                    gtab.append(('GO_%s' % gtype,string.join(ghtml,' ~ ')))
        if gtab:
            html += ['','<!-- ~~~~ %s GO tabs ~~~~ -->' % gene,tabberHTML('GO',gtab),
                     '<!-- ~~~~ End %s GO tabs ~~~~ -->' % gene,'']
    except: callobj.errorLog('seqDetailsHTML Error')
    ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    #print string.join(html,'\n')
    return string.join(html,'\n')
#########################################################################################################################
def checkHTML(hpage):   ### Checks for existence of complete HTML page
    '''Checks for existence of complete HTML page.'''
    if not os.path.exists(hpage): return False
    html = open(hpage,'r').read().upper()
    if html.find('<HTML') < 0 or html.find('</HTML>') < 0: return False
    return True
#########################################################################################################################
def stripTags(html,keeptags=[]):    ### Strips all HTML tag text from html code, except listed keeptags
    '''Strips all HTML tag text from html code, except listed keeptags.'''
    keeptags = string.split(string.join(keeptags).lower())
    tagsplit = string.split(html,'<')
    newhtml = tagsplit.pop(0)    
    while tagsplit:
        tagtxt = tagsplit.pop(0)
        tag = rje.matchExp('^\\\\?([A-Za-z0-9]+)',tagtxt)
        if tag and tag[0].lower() in keeptags: newhtml += '<%s' % tagtxt
        elif tagtxt.find('>') >= 0: newhtml += ' %s' % tagtxt[tagtxt.find('>')+1:]
    return string.replace(newhtml,'  ',' ')
#########################################################################################################################
def redirectToURL(url,refresh=0):    ### Returns HTML to redirect to URL. Should be within <HEAD>
    '''
    Returns HTML to redirect to URL.
    >> url:str = URL to redirect to.
    >> refresh:int [0] = Time to redirect
    '''
    return '<meta http-equiv="REFRESH" content="%d;url=%s">' % (refresh,url)
#########################################################################################################################
def tableToHTML(delimtext,delimit,tabwidth='100%',tdwidths=[],tdalign=[],valign='center',thead=True,border=1,tabid=''):    # Makes HTML Table
    '''
    Converts delimited plain text into an HTML table.
    >> delimtext:str = Delimited text to convert
    >> delimit:str = Text delimiter for conversion.
    >> tabwidth:str ['100%'] = width of table
    >> tdwidths:list [] = Optional list of widths of columns
    >> tdalign:list [] = Optional list of text alignment for columns
    >> valign:str ['center'] = Vertical text alignment for columns
    >> thead:bool [True] = Whether first row should use th rather than td
    >> border:int [1] = Table border strength
    >> tabid:str [''] = Table ID setting (for CSS formatting)
    '''
    ### [0] Setup Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if tabid: html = '<table width="%s" id="%s">\n' % (tabwidth,tabid)
    else: html = '<table width="%s" border=%d>\n' % (tabwidth,border)
    tablines = string.split(delimtext,'\n')
    ### [1] Header Row ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if thead:
        html += '<tr>\n'
        headtext = rje.readDelimit(tablines.pop(0),delimit)
        tw = tdwidths[0:]
        ta = tdalign[0:]
        while headtext:
            tag = 'th'
            if tw: tag += ' width="%s"' % tw.pop(0)
            if ta: tag += ' align=%s' % ta.pop(0)
            html += '<%s>%s</th>\n' % (tag,headtext.pop(0))
        html += '</tr>\n'
    ### [2] Main body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    while tablines:
        tabtext = rje.readDelimit(tablines.pop(0),delimit)
        if not tabtext: continue
        html += '<tr>\n'
        tw = tdwidths[0:]
        ta = tdalign[0:]
        while tabtext:
            tag = 'td valign=%s' % valign
            if tw: tag += ' width="%s"' % tw.pop(0)
            if ta: tag += ' align=%s' % ta.pop(0)
            html += '<%s>%s</td>\n' % (tag,tabtext.pop(0))
        html += '</tr>\n'
    ### [3] End table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    html += '</table>\n\n'
    return html
#########################################################################################################################
def dbTableToHTML(table,fields=[],datakeys=[],tabwidth='100%',tdwidths=[],tdalign=[],tdtitles={},valign='center',thead=True,border=1,tabid='',datasort={}):    # Makes HTML Table
    '''
    Converts delimited plain text into an HTML table.
    >> table:Database.Table object to convert
    >> fields:list [] = List of fields to output. Will use them all if empty.
    >> datakeys:list [] = List of entry data keys to output. Will use them all if empty.
    >> tabwidth:str ['100%'] = width of table
    >> tdwidths:list [] = Optional list of widths of columns
    >> tdalign:list [] = Optional list of text alignment for columns
    >> tdtitles:dict {} = Optional dictionary of {field:title text (for mouseover)}
    >> valign:str ['center'] = Vertical text alignment for columns
    >> thead:bool [True] = Whether first row should use th rather than td
    >> border:int [1] = Table border strength
    >> tabid:str [''] = Table ID setting (for CSS formatting)
    >> datasort:dict {'*':'string'} = Dictionary of field:type for stupidtable.js sorting. (* = default)
    '''
    try:### [0] Setup Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if tabid: html = '<table width="%s" id="%s">\n' % (tabwidth,tabid)
        else: html = '<table width="%s" border=%d>\n' % (tabwidth,border)
        if not fields: fields = table.fields()
        if not datakeys: datakeys = table.dataKeys()
        #table.debug(fields)
        ### [1] Header Row ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if thead:
            html += '<tr>\n'
            headtext = fields[0:]
            tw = tdwidths[0:]
            ta = tdalign[0:]
            while headtext:
                tag = 'th'
                if headtext[0] in datasort: tag += ' data-sort="%s"' % datasort[headtext[0]]
                elif '*' in datasort: tag += ' data-sort="%s"' % datasort['*']
                if tw: tag += ' width="%s"' % tw.pop(0)
                if ta: tag += ' align=%s' % ta.pop(0)
                if headtext[0] in tdtitles: tag += ' title="%s"' % tdtitles[headtext[0]]
                if headtext[0] in datasort or '*' in datasort: headtext[0] += ' <i class="fa fa-sort"></i>'
                html += '<%s>%s</th>\n' % (tag,headtext.pop(0))
            html += '</tr>\n'
        ### [2] Main body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for dkey in datakeys:
            entry = table.data(dkey)
            tabtext = []
            if not type(entry) == dict:
                table.warnLog('Table "%s" entry skipped for HTML - not a dictionary!: "%s"' % (dkey,entry))
                continue
            for field in fields: tabtext.append(entry[field])
            html += '<tr>\n'
            tw = tdwidths[0:]
            ta = tdalign[0:]
            while tabtext:
                tag = 'td valign=%s' % valign
                if tw: tag += ' width="%s"' % tw.pop(0)
                if ta: tag += ' align=%s' % ta.pop(0)
                html += '<%s>%s</td>\n' % (tag,tabtext.pop(0))
            html += '</tr>\n'
        ### [3] End table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        html += '</table>\n\n'
        return html
    except: table.errorLog('Problem with dbTableToHTML(%s)' % table.name())
#########################################################################################################################
def progStartHTML(callobj):     ### Returns <code><pre> HTML of start time, argument list and run directory
    '''Returns <code><pre> HTML of start time, argument list and run directory.'''
    ### [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    argcmd = rje.longCmd(sys.argv[1:])
    info = callobj.log.obj['Info']
    ### ~ [1] ~ Make and return HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    #html = '<code><pre>'
    html = '<pre><code>'
    html += 'Run Report for %s V%s: %s\n' % (info.program,info.version,time.asctime(time.localtime(info.start_time)))
    html += 'Run from directory: %s\n' % os.path.abspath(os.curdir)
    html += 'Commandline arguments: %s\n' % string.join(argcmd)
    #html += '</pre></code>\n\n'
    html += '</code></pre>\n\n'
    return html
#########################################################################################################################
### END OF SECTION III: MODULE METHODS                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: print 'Unexpected error during program setup:', sys.exc_info()[0]; return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: #HTML(mainlog,cmd_list).run()
        print rje_zen.Zen().wisdom(), '\n\n *** No standalone functionality! *** \n\n'


    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
