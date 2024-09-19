#!/usr/local/bin/python

# pic_html - HTML generator for picture websites
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
Module:       pic_html
Description:  HTML generator for picture websites
Version:      1.1
Last Edit:    23/11/05
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to take one or more folders of pictures and generate linked HTML pages for them.
    An accessory application for making thumbnails is needed. This should be used prior to running the program,
    or at the prompt.

    At present, this script will generate photos for an assumed hiearchical ordering. The current directory
    should be the one from which the program is run, and this will contain the pictures.htm file which links
    to the other pages. Within this there are then three hierarchies of folder:
    1. A 'Type' of picture group, e.g. Holidays or Animals etc.
    2. A 'Group' of pictures within the type, which may itself have multiple subfolders
    3. An 'Element' of the group of pictures, which consists of a single set of linked photos

    The 'Group' is the primary focus of this script as photos will generally be added for a particular event
    (a holiday for example). In addition to the folders for individual elements, each group should have a thumbnails
    folder which has all the thumbnails for that group. Each element has its own folder for which a set of linked
    web pages is generated. These are linked directly to the Pictures page.

    The pages for each element consist of two frames:
    1. A margin, for easy navigation to individual photos and thumbnails pages etc.
    2. The main frame where the photos are displayed, with their names and a short description.

    Each photo has its own page, with arrows cycling through the previous and next photos (in a circular fashion
    via the thumbnails page) with the title of the picture underneath and, if desired, a small description. Clicking
    on the picture will open it full size in another page. In addition to these, there is a thumbnails page, which
    links to the individual web-pages, and a downloads page which links direct to the photos themselves.

    If usehome=T and/or edithome=T are used (both are by default) the pictures home page will be read for existing
    descriptions and edited to contain links to all folders and descriptions in the picdesc file, i.e. the content
    from the two files plus any additional folders added, will be merged. When reading from the pictures home page,
    the script recognises Types, Groups and Elements of photos on a strict pattern recognition in the HTML:
        Type:   '^<HR><H2>(.+)</H2>'    => Type Description
        Group:  '^<P><B>(.+)</B>'       => Group Description
        Element:'^<A HREF="(\S+)/(\S+)/(\S+)/index.htm" TARGET="_top">(.+\S)</A>' => type,group,el, Element Description
    In each case, the number of photos may be extracted from the description with '\s+[(\d+)]$'. Be careful that no other
    lines will match these patterns. (The Group pattern may occur BEFORE the first Type pattern but NOT AFTER it.)

Commandline:
    ### File locations etc. ###
    pichome=FILE        : File name for pictures home page [pictures.htm]
    usehome=T/F         : Whether to extract descriptions etc. from the pictures home page [True]
    edithome=T/F        : Whether to regenerate the pictures home page [True]
    picdesc=FILE        : File with picture descriptions Type:Group:Element Description [pic_descriptions.txt]
    extlist=X,Y,..,Z    : List of acceptable picture file extensions [jpg,gif]
    picfolder=X         : Name of primary folder to look in for photos [*]
    thumbnails=X        : Folder containing thumbnails for all pictures of this type [thumbnails]
    thumbname=X         : Name to distinguish thumbnails from actual pictures (will rename) [_thumb]

    ### Webpage Appearance ###
    homeback=FILE   : File to use for picture home background [None]
    hometitle="X"   : Title to use for home web page ["Photo Page"]
    picback=FILE    : File to use for other page backgrounds [None]
    fontface=X      : Font to use for text [Comic Sans MS]
    larrow=FILE     : Image file for Left Arrow (pref GIF) [larrow.gif]
    rarrow=FILE     : Image file for Right Arrow (pref GIF) [rarrow.gif]
    picids=T/F      : Whether pictures have picture IDs 'ID - Name.*' [True]
    thumbheight=X   : Height of thumbnails (pixels) in Preview pages [120]

    ### Other ###
    clearhtml=T/F   : Delete existing HTML in picture folders (*.htm and *.html) - may overwrite anyway. [True]

Uses general modules: copy, glob, os, re, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy
import glob
import os
import re
import string
import sys
import time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#########################################################################################################################
### History
# 0.0 - Initial Compilation based on pic_html.plx.
# 1.0 - Full working version
# 1.1 - Improved Picture Group Ordering in thumbnail sites
#########################################################################################################################
### Major Functionality to Add
# [ ] : Decent Dealing of GIFs in thumbnails - makes an empty link at present!
# [ ] : Movies!
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'PIC_HTML'
        version = '1.1'
        last_edit = 'November 05'
        description = 'HTML generator for picture websites'
        author = 'Dr Richard J. Edwards.'
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print 'Problem making Info object.'
        raise
#########################################################################################################################
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
        print "Problem during initial setup."
        raise
#############################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
re_htmltype = '^<HR><H2>(.+)</H2>'  # => Type Description
re_htmlgroup = '^<P><B>(.+)</B>'    # => Group Description
re_htmlelement = '^<A HREF="(\S+)/(\S+)/(\S+)/index.htm" TARGET="_top">(.+\S)</A>'
    # => type,group,el, Element Description
re_htmlpicnum = '^(.*\S)\s+\[(\d+)\]$'    # => Number of photos from description
#########################################################################################################################
### SECTION II: CLASSES                                                                                                 #
#########################################################################################################################

#########################################################################################################################
### picHTML class:                                                                                                      #
#########################################################################################################################
class picHTML(rje.RJE_Object):     
    '''
    picHTML Class. Author: Rich Edwards (2005).

    Info:str
    - PicHome = File name for pictures home page [pictures.htm]
    - HomeTitle = Title to use for home web page ["Photo Page"]
    - PicDesc = File with picture descriptions Type:Group:Element Description [pic_descriptions.txt]
    - ExtList = List of acceptable picture file extensions [jpg,gif]
    - PicFolder = Name of primary folder to look in for photos [*]
    - ThumbNails = Folder containing thumbnails for all pictures of this type [thumbnails]
    - ThumbName = Name to distinguish thumbnails from actual pictures (will rename) [_thumb]
    - Background = File to use for page background [None]
    - HomeBack = File to use for picture home background [None]
    - FontFace = Font to use for text [Comic Sans MS]
    - LArrow = Image file for Left Arrow (pref GIF) [larrow.gif]
    - RArrow = Image file for Right Arrow (pref GIF) [rarrow.gif]
    
    Opt:boolean
    - UseHome  = Whether to extract descriptions etc. from the pictures home page [True]
    - EditHome = Whether to regenerate the pictures home page [True]
    - ClearHTML = Delete existing HTML in picture folders (*.htm and *.html) - may overwrite anyway. [False]
    - PicIDs = Whether pictures have picture IDs 'ID - Name.*' [True]

    Stat:numeric
    - ThumbHeight = Height of thumbnails (pixels) in Preview pages [120]

    Obj:RJE_Objects
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['PicHome','HomeTitle','PicDesc','ExtList','PicFolder','ThumbNails','ThumbName','HomeBack','Background','FontFace','LArrow','RArrow']
        - Opt:boolean ['ClearHTML','PicIDs','UseHome','EditHome']
        - Stats:float ['ThumbHeight']
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['PicHome','HomeTitle','PicDesc','ExtList','PicFolder','ThumbNails','ThumbName','Background','HomeBack','FontFace','LArrow','RArrow']
        self.optlist = ['ClearHTML','PicIDs','UseHome','EditHome']
        self.statlist = ['ThumbHeight']
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None)
        self.setInfo({'Name':'pic_html','PicHome':'pictures.htm','HomeTitle':'Photo Page',
                      'PicDesc':'pic_descriptions.txt','ExtList':'jpg,gif','PicFolder':'*',
                      'ThumbNails':'thumbnails','ThumbName':'_thumb'})
        self.setInfo({'FontFace':'Comic Sans MS','LArrow':'larrow.gif','RArrow':'rarrow.gif'})
        self.opt['ClearHTML'] = False
        self.stat['ThumbHeight'] = 120
        ### <c> ### Other Attributes
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for i in range(len(self.cmd_list)):
            cmd = self.cmd_list[i]
            try:
                ### <a> ### General Options
                self._generalCmd(cmd)
                ### <b> ### Class Options
                self._cmdRead(cmd,type='info',att='PicHome')  
                self._cmdRead(cmd,type='info',att='PicDesc')  
                self._cmdRead(cmd,type='info',att='ExtList')  
                self._cmdRead(cmd,type='info',att='PicFolder')  
                self._cmdRead(cmd,type='info',att='ThumbNails')
                self._cmdRead(cmd,type='info',att='ThumbName') 
                self._cmdRead(cmd,type='info',att='Background',arg='picback')  
                self._cmdRead(cmd,type='info',att='Background')  
                self._cmdRead(cmd,type='info',att='HomeBack')  
                self._cmdRead(cmd,type='info',att='HomeTitle')  
                self._cmdRead(cmd,type='info',att='FontFace')  
                self._cmdRead(cmd,type='info',att='LArrow')  
                self._cmdRead(cmd,type='info',att='RArrow')  
                self._cmdRead(cmd,type='opt',att='UseHome')
                self._cmdRead(cmd,type='opt',att='EditHome')
                self._cmdRead(cmd,type='opt',att='ClearHTML')
                self._cmdRead(cmd,type='opt',att='PicIDs')  
                self._cmdRead(cmd,type='stat',att='ThumbHeight')
                ### Special
                if cmd.find('hometitle="') == 0:
                    x = i + 0
                    hometitle = [cmd[len('hometitle="'):]]
                    while self.cmd_list[x][-1] != '"' and x < (len(self.cmd_list) - 1):
                        x += 1
                        hometitle.append(self.cmd_list[x])
                    self.info['HomeTitle'] = rje.join(hometitle,' ')[:-1]
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Methods                                                                                            #
#########################################################################################################################
    def run(self):      ### Main run method
        '''
        Generic method. Add description here (and arguments.)
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            ## <a> ## User choice
            self.verbose(0,2,'\n\n%s' % self.details(),0)
            if self.stat['Interactive'] >= 0 and not rje.yesNo('Proceed with these settings?'):
                self.edit()
            extlist = rje.split(self.info['ExtList'],',')
            ## <b> ## Types, Groups and Elements
            folderpath = rje.split(rje.makePath(self.info['PicFolder'],wholepath=True),os.sep)
            if len(folderpath) > 3:
                self.verbose(0,1,'Cannot process beyond Elements (*/*/* at this stage!',1)
                return
            # Type 
            typelist = []   # List of Type folders to process
            groupdic = {}   # List of Group folders to process for each Type key
            elementdic  = {} # List of Element folders to process for each Type:Group key
            pixdic = {}     # List of pictures in Type:Group:Element
            grouppath = '*'
            elpath = '*'
            if len(folderpath) == 3:
                typelist = folderpath[0:1]
                [grouppath,elpath] = folderpath[1:3]
            elif len(folderpath) == 2:
                typelist = folderpath[0:1]
                grouppath = folderpath[1]
            else:
                typelist = rje.subDir(folderpath[0],[self.info['ThumbNails']])
            # Groups & Elements
            for type in typelist:
                groupdic[type] = rje.subDir(rje.join([type,grouppath],os.sep),[self.info['ThumbNails']])
                elementdic[type] = {}
                pixdic[type] = {}
                for group in groupdic[type]:
                    if self.opt['ClearHTML']:
                        html_die = glob.glob(rje.join([type,group,'*.htm'],os.sep)) + glob.glob(rje.join([type,group,'*.html'],os.sep))
                        for file in html_die:
                            os.unlink(file)
                    elementdic[type][group] = rje.subDir(rje.join([type,group,elpath],os.sep),[self.info['ThumbNails']])
                    pixdic[type][group] = {}
                    # Pictures
                    for el in elementdic[type][group]:
                        if self.opt['ClearHTML']:
                            html_die = glob.glob(rje.join([type,group,el,'*.htm'],os.sep)) + glob.glob(rje.join([type,group,el,'*.html'],os.sep))
                            for file in html_die:
                                os.unlink(file)
                        pixdic[type][group][el] = []
                        for ext in extlist:
                            piclist = glob.glob(rje.join([type,group,el,'*.%s' % ext],os.sep))
                            for pic in piclist:
                                pixdic[type][group][el].append(rje.split(pic,os.sep)[-1])
            self.sumPix(pixdic) #!# Add option to accept/reject folders
            self.deBug(pixdic)
            ## <c> ## Descriptions
            eldesc = self._makeDesc(pixdic) # Dictionary of element descriptions

            ### <1> ### Thumbnails 
            for type in rje.sortKeys(pixdic):
                self.verbose(0,4,'\nChecking thumbnails for %s: ' % type,0)
                for group in rje.sortKeys(pixdic[type]):
                    while 1:
                        tdir = rje.join([type,group,self.info['ThumbNails']],os.sep)
                        renamex = 0
                        missing = []
                        if not rje.checkForFile(tdir):
                            os.mkdir(tdir)
                            self.verbose(0,2,'Directory %s Made:...' % tdir,0)
                            for el in pixdic[type][group]:
                                missing += pixdic[type][group][el][0:]
                        else:
                            for el in pixdic[type][group]:
                                for pix in pixdic[type][group][el]:
                                    pic = rje.join([tdir,pix],os.sep)
                                    tpic = self._thumbName(pic)
                                    if rje.checkForFile(pic):
                                        if rje.checkForFile(tpic):
                                            self.log.errorLog('Cannot rename %s as %s already exists! Check that IDs are unique!',False,True)
                                        else:                                            
                                            os.rename(pic,tpic)
                                            renamex += 1
                                    if not rje.checkForFile(tpic):  # Missing!
                                        if tpic[-4:].lower() == '.gif':
                                            self.verbose(0,1,'\n!!! GIF thumbnails not supported. Make thumbnail for %s as %s !!!' % (rje.join([type,group,el,pix],os.sep),tpic),2)
                                        else:
                                            missing.append(pix)
                        self.verbose(0,2,'Group %s: %d renamed, %d missing...' % (group,renamex,len(missing)),0)                            
                        if len(missing) > 0:    # Missing thumbnail pictures
                            self.verbose(0,4,'\n\nWARNING! The following %s thumbnails are missing:\n\n%s' % (group,missing),2)
                            if rje.yesNo('Generate thumbnails and recheck (y) or skip with broken links (n)?:'):
                                rje.choice('Generate thumbnails in %s and hit ENTER when ready.' % tdir)
                            else:
                                break
                        else:
                            break
                
            ### <2> ### Create Linked webpages for each element
            for type in rje.sortKeys(pixdic):
                self.verbose(0,4,'\nMaking webpages for %s: ' % type,0)
                for group in rje.sortKeys(pixdic[type]):
                    ## Group index.htm (Frames) ##
                    HTML = open(rje.join([type,group,'index.htm'],os.sep),'w')
                    HTML.write('<HTML><HEAD><TITLE>%s</TITLE></HEAD>\n\n' % eldesc[type][group]['Desc'])
                    HTML.write('<FRAMESET COLS=25%%,75%%>\n<FRAME SRC="%s_margin.htm" NAME="Margin">\n' % group)
                    HTML.write('<FRAME SRC="%s_thumbnails.htm" NAME="Main">\n' % group)
                    HTML.write("<NOFRAMES>\n<BODY>\n<H1>You don't support frames!</H1>\n")
                    HTML.write('<P>You can enter the site without frames <A HREF="%s_thumbnails.htm">here</A> ' % group)
                    HTML.write('but the site may not work properly. Sorry about that. <I>Go get a browser that supports frames!</I></P>\n')
                    HTML.write('</BODY>\n</NOFRAMES>\n</FRAMESET>\n</HTML>\n')
                    HTML.close()
                    self.verbose(2,3,'Index made...',0)

                    ## Group margin.htm
                    HTML = open(rje.join([type,group,'%s_margin.htm' % group],os.sep),'w')
                    HTML.write('<HTML><HEAD><TITLE>%s Margin</TITLE></HEAD>\n\n' % eldesc[type][group]['Desc'])
                    if self.info['Background'] != 'None':
                        HTML.write('<BODY BACKGROUND="../../%s" WIDTH="100%%">\n\n' % self.info['Background'])
                    HTML.write('<FONT FACE="%s" SIZE=-1><CENTER>\n' % self.info['FontFace'])
                    HTML.write('<P><A HREF="../../%s" TARGET="_top">Picture Home</A><HR>\n' % self.info['PicHome'])
                    HTML.write('<A HREF="%s_thumbnails.htm" TARGET="Main">Thumbnails</A><HR>\n' % group)
                    HTML.write('</CENTER><UL>')
                    for el in rje.sortKeys(pixdic[type][group]):
                        HTML.write('<LI><A HREF="%s/index.htm" TARGET="_top">%s</A><BR>\n' % (el,eldesc[type][group][el]))
                    HTML.write('</UL></FONT></BODY>\n</HTML>\n')
                    HTML.close()
                    self.verbose(2,3,'Margin made...',0)

                    ## Group thumbnails.htm
                    HTML = open(rje.join([type,group,'%s_thumbnails.htm' % group],os.sep),'w')
                    HTML.write('<HTML><HEAD><TITLE>%s Thumbnails</TITLE></HEAD>\n\n' % eldesc[type][group]['Desc'])
                    if self.info['Background'] != 'None':
                        HTML.write('<BODY BACKGROUND="../../%s" WIDTH="100%%">\n\n' % self.info['Background'])
                    HTML.write('<H2><CENTER>%s</H2>\n' % eldesc[type][group]['Desc'])
                    HTML.write('<FONT FACE="%s" SIZE=-1>\n' % self.info['FontFace'])
                    for el in rje.sortKeys(pixdic[type][group]):
                        HTML.write('<HR><H3><A HREF="%s/index.htm" TARGET="_top">%s</A></H3>' % (el,eldesc[type][group][el]))
                        for pic in pixdic[type][group][el]:
                            tpic = self._thumbName(pic)
                            HTML.write('<A HREF="%s/%s" TARGET="_Main">' % (el,pic))
                            HTML.write('<IMG SRC="%s/%s" ALT="%s" ' % (self.info['ThumbNails'],tpic,self._picName(pic)))
                            HTML.write('BORDER=2 HEIGHT=%d></A>\n' % self.stat['ThumbHeight'])
                    HTML.write('</FONT></BODY>\n</HTML>\n')
                    HTML.close()
                    self.verbose(2,3,'Preview made...',0)

                    self.verbose(0,4,'%s...' % group,0)
                    for el in rje.sortKeys(pixdic[type][group]):
                        ## Element index.htm (Frames) ##
                        HTML = open(rje.join([type,group,el,'index.htm'],os.sep),'w')
                        HTML.write('<HTML><HEAD><TITLE>%s</TITLE></HEAD>\n\n' % eldesc[type][group][el])
                        HTML.write('<FRAMESET COLS=25%%,75%%>\n<FRAME SRC="%s_margin.htm" NAME="Margin">\n' % el)
                        HTML.write('<FRAME SRC="%s_thumbnails.htm" NAME="Main">\n' % el)
                        HTML.write("<NOFRAMES>\n<BODY>\n<H1>You don't support frames!</H1>\n")
                        HTML.write('<P>You can enter the site without frames <A HREF="%s_thumbnails.htm">here</A> ' % el)
                        HTML.write('but the site may not work properly. Sorry about that. <I>Go get a browser that supports frames!</I></P>\n')
                        HTML.write('</BODY>\n</NOFRAMES>\n</FRAMESET>\n</HTML>\n')
                        HTML.close()
                        self.verbose(2,3,'Index made...',0)

                        ## Element margin.htm
                        HTML = open(rje.join([type,group,el,'%s_margin.htm' % el],os.sep),'w')
                        HTML.write('<HTML><HEAD><TITLE>%s Margin</TITLE></HEAD>\n\n' % eldesc[type][group][el])
                        if self.info['Background'] != 'None':
                            HTML.write('<BODY BACKGROUND="../../../%s" WIDTH="100%%">\n\n' % self.info['Background'])
                        HTML.write('<FONT FACE="%s" SIZE=-1><CENTER>\n' % self.info['FontFace'])
                        HTML.write('<P><A HREF="../../../%s" TARGET="_top">Picture Home</A><HR>\n' % self.info['PicHome'])
                        HTML.write('<A HREF="../index.htm" TARGET="_top">All %s</A><HR>\n' % eldesc[type][group]['Desc'])
                        HTML.write('<A HREF="%s_downloads.htm" TARGET="Main">Downloads</A><HR>\n' % el)
                        HTML.write('<A HREF="%s_thumbnails.htm" TARGET="Main">Thumbnails</A><HR>\n' % el)
                        HTML.write('</CENTER><UL>')
                        for pic in pixdic[type][group][el]:
                            HTML.write('<LI><A HREF="%s.htm" TARGET="Main">%s</A><BR>\n' % (self._picID(pic),self._picName(pic)))
                        HTML.write('</UL></FONT></BODY>\n</HTML>\n')
                        HTML.close()
                        self.verbose(2,3,'Margin made...',0)

                        ## Element thumbnails.htm
                        HTML = open(rje.join([type,group,el,'%s_thumbnails.htm' % el],os.sep),'w')
                        HTML.write('<HTML><HEAD><TITLE>%s Thumbnails</TITLE></HEAD>\n\n' % eldesc[type][group][el])
                        if self.info['Background'] != 'None':
                            HTML.write('<BODY BACKGROUND="../../../%s" WIDTH="100%%">\n\n' % self.info['Background'])
                        HTML.write('<H2><CENTER>%s (Thumbnails)</H2>\n' % eldesc[type][group][el])
                        HTML.write('<FONT FACE="%s" SIZE=-1>\n' % self.info['FontFace'])
                        for pic in pixdic[type][group][el]:
                            tpic = self._thumbName(pic)
                            HTML.write('<A HREF="%s.htm">' % self._picID(pic))
                            HTML.write('<IMG SRC="../%s/%s" ALT="%s" ' % (self.info['ThumbNails'],tpic,self._picName(pic)))
                            HTML.write('BORDER=3 HEIGHT=%d></A>\n' % self.stat['ThumbHeight'])
                        HTML.write('</FONT></BODY>\n</HTML>\n')
                        HTML.close()
                        self.verbose(2,3,'Preview made...',0)

                        ## Element downloads.htm
                        HTML = open(rje.join([type,group,el,'%s_downloads.htm' % el],os.sep),'w')
                        HTML.write('<HTML><HEAD><TITLE>%s Downloads</TITLE></HEAD>\n\n' % eldesc[type][group][el])
                        if self.info['Background'] != 'None':
                            HTML.write('<BODY BACKGROUND="../../../%s" WIDTH="100%%">\n\n' % self.info['Background'])
                        HTML.write('<H2><CENTER>%s (Downloads)</H2>\n' % eldesc[type][group][el])
                        HTML.write('<FONT FACE="%s" SIZE=-1>\n' % self.info['FontFace'])
                        HTML.write('<P><I>This page links directly to the image files, which will be opened in new windows.<BR>\n')
                        HTML.write('Click <A HREF="%s_full_images.htm">here</A> to open all %d full-size images on one page.</I></P>\n\n' % (el,len(pixdic[type][group][el])))
                        for pic in pixdic[type][group][el]:
                            tpic = self._thumbName(pic)
                            HTML.write('<A HREF="%s" TARGET="_blank">' % pic)
                            HTML.write('<IMG SRC="../%s/%s" ALT="%s" ' % (self.info['ThumbNails'],tpic,pic))
                            HTML.write('BORDER=3 HEIGHT=%d></A>\n' % self.stat['ThumbHeight'])
                        HTML.write('</FONT></BODY>\n</HTML>\n')
                        HTML.close()
                        self.verbose(2,3,'Downloads made...',0)

                        ## Element downloads.htm
                        HTML = open(rje.join([type,group,el,'%s_full_images.htm' % el],os.sep),'w');
                        HTML.write('<HTML><HEAD><TITLE>%s Full Images</TITLE></HEAD>\n\n' % eldesc[type][group][el])
                        if self.info['Background'] != 'None':
                            HTML.write('<BODY BACKGROUND="../../../%s" WIDTH="100%%">\n\n' % self.info['Background'])
                        HTML.write('<H2><CENTER>%s (Full Size)</H2>\n' % eldesc[type][group][el])
                        HTML.write('<FONT FACE="%s" SIZE=-1>\n' % self.info['FontFace'])
                        HTML.write('<P><I>The images on this page are <B>full size files</B>, <U>not</U> thumbnails.</I></P>\n\n')
                        for pic in pixdic[type][group][el]:
                            HTML.write('<A HREF="%s"  TARGET="_blank">' % pic)
                            HTML.write('<IMG SRC="%s" ALT="%s" ' % (pic,pic))
                            HTML.write('BORDER=3 HEIGHT=%d></A>\n' % int(1.5 * self.stat['ThumbHeight']))
                        HTML.write('</FONT></BODY>\n</HTML>\n')
                        HTML.close()
                        self.verbose(2,3,'Full_Image made...',0)

                        ### Individual Element Pages ###
                        for pic in pixdic[type][group][el]:
                            HTML = open(rje.join([type,group,el,'%s.htm' % self._picID(pic)],os.sep),'w')
                            HTML.write('<HTML><HEAD><TITLE>%s</TITLE></HEAD>\n\n' % self._picName(pic))
                            if self.info['Background'] != 'None':
                                HTML.write('<BODY BACKGROUND="../../../%s" WIDTH="100%%">\n\n' % self.info['Background'])
                            HTML.write('<FONT FACE="%s" SIZE=-1><CENTER>\n' % self.info['FontFace'])
                            HTML.write('<IMG SRC="%s" ALT="%s" BORDER=0 HEIGHT=85%%><HR>' % (pic,pic))
                            i = pixdic[type][group][el].index(pic)
                            HTML.write('<TABLE WIDTH=100%><TR><TD WIDTH=15% ALIGN=CENTER>\n')
                            if i == 0:
                                #HTML.write('<A HREF="%s_thumbnails.htm"><IMG SRC="../../../%s" ALT="Thumbnails" BORDER=0 WIDTH=50></A>\n' % (el,self.info['LArrow']))
                                HTML.write('<A HREF="%s_thumbnails.htm">Preview</A>\n' % el)
                            else:
                                HTML.write('<A HREF="%s.htm"><IMG SRC="../../../%s" ALT="%s" BORDER=0 WIDTH=50></A>\n' % (self._picID(pixdic[type][group][el][i-1]),self.info['LArrow'],self._picName(pixdic[type][group][el][i-1])))
                            HTML.write('</TD><TD WIDTH=70%% ALIGN=CENTER>%s\n' % self._picName(pic))
                            HTML.write('</TD><TD WIDTH=15% ALIGN=CENTER>\n')
                            if i == len(pixdic[type][group][el]) - 1:
                                #HTML.write('<A HREF="%s_thumbnails.htm"><IMG SRC="../../../%s" ALT="Thumbnails" BORDER=0 WIDTH=50></A>\n' % (el,self.info['RArrow']))
                                HTML.write('<A HREF="%s_thumbnails.htm">Preview</A>\n' % el)
                            else:
                                HTML.write('<A HREF="%s.htm"><IMG SRC="../../../%s" ALT="%s" BORDER=0 WIDTH=50></A>\n' % (self._picID(pixdic[type][group][el][i+1]),self.info['RArrow'],self._picName(pixdic[type][group][el][i+1])))
                            HTML.write('</TD></TR></TABLE>\n</FONT></BODY></HTML>\n')
                            HTML.close()
                            self.verbose(2,3,'%s...' % pic,0)

                self.verbose(0,1,'Done!',1)
            
            ### <3> ### Update Main Page
            if not self.opt['EditHome']:
                return
            self.verbose(0,4,'\nReMaking picture home webpage (%s): ' % self.info['PicHome'],0)
            ## Heading etc. ##
            HTML = open(self.info['PicHome'],'w')
            HTML.write('<HTML>\n<HEAD><TITLE>%s</TITLE></HEAD>\n\n' % self.info['HomeTitle'])
            if self.info['HomeBack'] != 'None':
                HTML.write('<BODY BACKGROUND="%s" WIDTH="100%%">\n\n' % self.info['HomeBack'])
            HTML.write('<FONT FACE="%s">\n' % self.info['FontFace'])            
            ## Description ##
            HTML.write('<H2><CENTER>%s</CENTER></H2>\n' % self.info['HomeTitle'])
            HTML.write('''<P>The photos on this page are broken down into basic categories, some of which have
            subcategories and stuff. Each 'page' of photos consists of a main frame in which the pictures appear
            and a link strip down the side for quick navigation to specific pics. (<B>NB.</B> The size of the
            photo can be increased/decreased by making the side strip thinner/wider. If saved to disk, the picture
            will probably be a different size to that viewed in this website - usually bigger.)</P>
            <P>Alternative navigation is provided by arrows below each picture to advance/retreat one picture at
            a time or by clicking on pictures on the 'Thumbnails' page. The 'Downloads page links directly to the
            pictures themselves (rather than to HTML) for easy downloading.</P>
            <P>(<B>NB.</B> Thumbnails are not currently available for GIFs and Movie files. Movie file webpages are
            also somewhat dodgy - go to the 'Downloads' page and click on the links to view.)</P>
            ''')
            ### Pictures
            for type in rje.sortKeys(eldesc):
                HTML.write('\n<HR><H2>%s</H2>\n' % eldesc[type]['Group']['Desc'])
                for group in rje.sortKeys(eldesc[type]):
                    if group in ['Group','Num']:
                        continue
                    HTML.write('<P><B>%s</B> [<A HREF="%s/%s/index.htm">All</A>]<BR>\n' % (eldesc[type][group]['Desc'],type,group))
                    within_el = False
                    for el in rje.sortKeys(eldesc[type][group]):
                        if el == 'Desc':
                            continue
                        if within_el:
                            HTML.write(' ~ ')
                        if type in pixdic.keys() and group in pixdic[type].keys() and el in pixdic[type][group].keys() and '%d' % len(pixdic[type][group][el]) != self._descNum(eldesc[type][group][el]): # New!
                            self.log.printLog('#NEW','%s/%s/%s => [%d] pics vs %s: New!' % (type,group,el,len(pixdic[type][group][el]),eldesc[type][group][el]))
                            HTML.write('<FONT COLOR=RED><I>New! </I></FONT>\n')
                            eldesc[type][group][el] = '%s [%d]' % (self._pureDesc(eldesc[type][group][el]),len(pixdic[type][group][el]))
                        HTML.write('<A HREF="%s/%s/%s/index.htm" TARGET="_top">%s</A>\n' % (type,group,el,eldesc[type][group][el]))
                        within_el = True
                    HTML.write('</P>\n\n')

            ## Update ##
            HTML.write('\n<HR><FONT COLOR=BLUE SIZE=2><B>Last modified:</B> %s<BR>\n' % time.asctime(time.localtime(time.time())))
            HTML.write('This site was generated automatically by the pic_html Python script (R Edwards 2005).</FONT></P>\n\n')
            ## End ##
            HTML.write('</FONT></BODY></HTML>\n')
            HTML.close()
            self.verbose(0,0,'Done!',1)

        except:
            self.log.errorLog('Error in _method(%s)' % _stage,printerror=True,quitchoice=True)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def _thumbName(self,pic):   ### Returns thumbnail name for pic
        '''
        Returns thumbnail name for pic.
        pic:str = picture name
        '''
        tpic = pic
        if self.opt['PicIDs'] and re.search('^(\S+) - (\S.+)$',pic):
            tpic = rje.matchExp('^(\S+) - (\S.+)$',pic)[0] + pic[-4:]
        tpic = tpic[:-4] + self.info['ThumbName'] + tpic[-4:]
        return tpic
#########################################################################################################################
    def sumPix(self,pixdic):      ### Prints summary of chosen elements
        '''
        Prints summary of chosen picture folders and contents.
        >> pixdic:List of pictures in Type:Group:Element
        '''
        try:
            (px,ex,gx,tx) = (0,0,0,0)
            self.verbose(0,4,'',1)
            for type in rje.sortKeys(pixdic):
                tx += 1
                self.verbose(0,4,'\n%s' % type,1)
                for group in rje.sortKeys(pixdic[type]):
                    gx += 1
                    self.verbose(0,4,'|-%s' % group,1)
                    for el in rje.sortKeys(pixdic[type][group]):
                        ex += 1
                        px += len(pixdic[type][group][el])
                        self.verbose(0,4,'  |-%s (%d pics)' % (el,len(pixdic[type][group][el])),1)
            self.verbose(0,0,'\n%d Pictures from %d Elements in %d Groups across %d Types.' % (px,ex,gx,tx),2)
            #!# Add option to accept/reject folders if rje.yesNo('Use all folders?):
        except:
            self.log.errorLog('Error in sumElements()',printerror=True,quitchoice=False)
#########################################################################################################################
    def _getDesc(self,folder,desc=None):  ### Returns user-input decription for folder
        '''
        Returns user-input decription for folder.
        >> folder:str = Name of folder
        >> desc:str = Default description [folder if None]
        << desc:str = description returned
        '''
        try:
            if not desc:
                desc = rje.split(folder,'/')[-1]
                desc = '%s%s' % (desc[0].upper(),desc[1:].lower())
            default = ''
            while desc != default and self.stat['Interactive'] > 0:
                default = desc
                desc = rje.choice('Description for folder "%s":' % folder,default)
            return desc
        except:
            self.log.errorLog('Error in _getDesc()',printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def _makeDesc(self,pixdic):  ### Makes Dictionary of element descriptions
        '''
        Makes Dictionary of element descriptions.
        >> pixdic:List of pictures in Type:Group:Element
        '''
        try:
            ### Load previous descriptions ###
            descdic = {}
            ## Load from text file ##
            if rje.checkForFile(self.info['PicDesc']):
                DESC = open(self.info['PicDesc'],'r')
                desclines = DESC.readlines()
                DESC.close()
                for line in desclines:
                    if re.search('^(\S+)\s+(\S.+)$',line):
                        (keys,desc) = rje.matchExp('^(\S+)\s+(\S.+)$',line)
                        keys = rje.split(keys,':')
                        if keys[-1] == '':
                            keys = keys[:-1]
                        if len(keys) != 3:
                            self.log.errorLog('Cannot get description for %s. Must be Type:Group:Element' % keys,False,False)
                            continue
                        if keys[0] not in descdic.keys():
                            descdic[keys[0]] = {keys[1]:{keys[2]:desc}}
                        elif keys[1] not in descdic[keys[0]].keys():
                            descdic[keys[0]][keys[1]] = {keys[2]:desc}
                        else:
                            descdic[keys[0]][keys[1]][keys[2]] = desc
            ## Load from pictures home file ##
            #print self.opt['UseHome'], rje.checkForFile(self.info['PicHome'])
            if self.opt['UseHome'] and rje.checkForFile(self.info['PicHome']):
                HOME = open(self.info['PicHome'],'r')
                homelines = HOME.readlines()
                HOME.close()
                (type,group,el,typedesc,groupdesc,eldesc) = ['None'] * 6
                for line in homelines:
                    #print line, rje.matchExp(re_htmltype,line), rje.matchExp(re_htmlgroup,line), rje.matchExp(re_htmlelement,line)
                    if re.search(re_htmltype,line): # New Type
                        typedesc = rje.matchExp(re_htmltype,line)[0]
                        groupdesc = 'None'
                    elif re.search(re_htmlgroup,line) and typedesc != 'None': # New Group
                        groupdesc = rje.matchExp(re_htmlgroup,line)[0]
                    elif re.search(re_htmlelement,line) and groupdesc != 'None':    # New Element
                        (type,group,el,eldesc) = rje.matchExp(re_htmlelement,line)
                        if type not in descdic.keys():
                            descdic[type] = {'Group':{'Desc':typedesc},group:{'Desc':groupdesc,el:eldesc}}
                        elif group not in descdic[type].keys():
                            descdic[type]['Group'] = {'Desc':typedesc}
                            descdic[type][group] = {'Desc':groupdesc,el:eldesc}
                        else:
                            descdic[type][group][el] = eldesc

            ### Get/Update relevant Descriptions
            for type in rje.sortKeys(pixdic):
                ## Type ##
                if type not in descdic.keys():
                    descdic[type] = {}
                if 'Group' not in descdic[type].keys():
                    descdic[type]['Group'] = {}
                if 'Desc' not in descdic[type]['Group'].keys():
                    descdic[type]['Group']['Desc'] = None
                #print type, descdic[type]['Group']['Desc']
                descdic[type]['Group']['Desc'] = self._getDesc(type,descdic[type]['Group']['Desc'])
                ## Group ##
                for group in rje.sortKeys(pixdic[type]):
                    if group not in descdic[type].keys():
                        descdic[type][group] = {}
                    if 'Desc' not in descdic[type][group].keys():
                        descdic[type][group]['Desc'] = None
                    #print type, group, descdic[type][group]['Desc']
                    descdic[type][group]['Desc'] = self._getDesc('%s/%s' % (type,group),descdic[type][group]['Desc'])
                ## Element ##
                    for el in rje.sortKeys(pixdic[type][group]):
                        if el not in descdic[type][group].keys():
                            descdic[type][group][el] = None
                        #print type, group, el, descdic[type][group][el]
                        descdic[type][group][el] = self._getDesc('%s/%s/%s' % (type,group,el),descdic[type][group][el])

            ### Save updated list ###
            DESC = open(self.info['PicDesc'],'w')
            for type in rje.sortKeys(descdic):
                for group in rje.sortKeys(descdic[type]):
                    for el in rje.sortKeys(descdic[type][group]):
                        DESC.write('%s:%s:%s: %s\n' % (type,group,el,self._pureDesc(descdic[type][group][el])))
            DESC.close()
            return descdic            
        except:
            self.log.errorLog('Error in _makeDesc()',printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def _pureDesc(self,desc):   ### Returns pure description minus any numbers of pictures
        '''
        Returns pure description minus any numbers of pictures.
        >> desc:str = Description with possible [number]
        '''
        if re.search(re_htmlpicnum,desc):
            return rje.matchExp(re_htmlpicnum,desc)[0]
        return desc
#########################################################################################################################
    def _descNum(self,desc):   ### Returns any numbers of pictures in description as string
        '''
        Returns any numbers of pictures in description as string.
        >> desc:str = Description with possible [number]
        '''
        if re.search(re_htmlpicnum,desc):
            return rje.matchExp(re_htmlpicnum,desc)[1]
        return ''
#########################################################################################################################
    def _picName(self,pic): ### Returns picture name from filename
        '''
        Returns picture name from filename.
        >> pic:str = File name
        '''
        picname = pic[:-4]  # Trim extension
        if self.opt['PicIDs'] and re.search('^(\S+) - (\S.+)$',picname):
            match = rje.matchExp('^(\S+) - (\S.+)$',picname)
            picname = '%s (%s)' % (match[1],match[0])
        return picname
#########################################################################################################################
    def _picID(self,pic): ### Returns picture ID from filename
        '''
        Returns picture ID from filename.
        >> pic:str = File name
        '''
        if self.opt['PicIDs'] and re.search('^(\S+) - (\S.+)$',pic):
            return rje.matchExp('^(\S+) - (\S.+)$',pic)[0]
        return pic
#########################################################################################################################
### End of picHTML                                                                                                      #
#########################################################################################################################
       
#########################################################################################################################
## End of SECTION II                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    try:
        ### <0>  ### Basic Setup of Program
        [info,out,mainlog,cmd_list] = setupProgram()
        
        ### <1> ### Rest...
        pix = picHTML(mainlog,cmd_list)
        pix.run()

        ### <X> ### End
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        mainlog.errorLog("User terminated.\n")
    except:
        print "Unexpected error:", sys.exc_info()[0]
    mainlog.printLog('#LOG', "%s V:%s End: %s\n" % (info.program, info.version, time.asctime(time.localtime(time.time()))), 1)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except SystemExit:
        sys.exit()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
