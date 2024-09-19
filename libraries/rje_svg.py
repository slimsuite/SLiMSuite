#!/usr/bin/python

# See below for name and description
# Copyright (C) 2009 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_svg
Description:  RJE SVG Module
Version:      0.0
Last Edit:    31/12/10
Copyright (C) 2009  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    ### ~~~ INPUT ~~~ ###
    col=LIST    : Replace standard colour listing (mixed Hex and RGB) []

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_zen
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
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : List here
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_SVG', '0.0', 'December 2010', '2010')
    description = 'RJE SVG Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print('\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print('Major Problem with cmdHelp()')
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
    except: print('Problem during initial setup.'); raise
#########################################################################################################################
### SVG Colours ###
soton = ["#014359",	# Marine Blue
         "#007275",	# Turquoise
         "#0A96A9",	# Light blue
         "#323D43",	# Dark grey
         "#979E45",	# Gold
         "#BBBBBB","#9BA3A6",	# Metallic ~ 1y End ~
         "#653A28",	# Dark Red
         "#531F44",	# Maroon
         "#A67891",	# Dusky Pink
         "#B2699F",	# Pink
         "#CCDAEA",	# Pale light blue ~ 2y1 End ~
         "#8A412B",	# Brown
         "#AB1210",	# Brick red
         "#F00F2C",	# Red
         "#FE3E14",	# Orange
         "#FFB300",	# Yellow ~ 2y2 End ~
         "#4F5A20",	# Green
         "#91BA91",	# Pale green
         "#BDB68A",	# Pale brown
         "#8F9E94"]	# Green/grey ~ End   
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class SVG(rje.RJE_Object):     
    '''
    SVG Class. Author: Rich Edwards (2010).

    Info:str
    
    Opt:boolean

    Stat:numeric

    List:list
    Col = Standard colour listing (mixed Hex and RGB)
    
    Dict:dictionary
    Col = Custom colour dictionaries. 

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = []
        self.optlist = []
        self.statlist = []
        self.listlist = ['Col']
        self.dictlist = ['Col']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdReadList(cmd,'list',['Col'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setupCol()
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def cwCol(self,aa,consaa):
        '''Returns ClustalW colour given aa and ClustalW conservation aa grouping.'''
        '''        
            if(sum(aa == c("X","x",".")) > 0){ return(soton$col[6]) }
            if(aa == "-"){ return("black") }
            if(aa == "G"){ return(rgbindex$ORANGE) }
            if(aa == "P"){ return(rgbindex$YELLOW) }
            if(aa == "T" & cwcontains(c("t","T","S","%","#"),consaa)){ return(rgbindex$GREEN) }
            if(aa == "S" & cwcontains(c("t","T","S","#"),consaa)){ return(rgbindex$GREEN) }
            if(aa == "N" & cwcontains(c("n","N","D"),consaa)){ return(rgbindex$GREEN) }
            if(aa == "Q" & cwcontains(c("q","Q","E","K","R","+"),consaa)){ return(rgbindex$GREEN) }
            if(sum(aa == c("W","L","V","I","M","F")) > 0 & cwcontains(c("%","#","A","C","F","H","I","L","M","V","W","Y","P","p"),consaa)){ return(rgbindex$BLUE) }
            if(aa == "A" & cwcontains(c("%","#","A","C","F","H","I","L","M","V","W","Y","P","p","T","S","s","G"),consaa)){ return(rgbindex$BLUE) }
            if(aa == "C" & cwcontains(c("%","#","A","S","F","H","I","L","M","V","W","Y","P","p"),consaa)){ return(rgbindex$BLUE) }
            if(aa == "C" & cwcontains(c("C"),consaa)){ return(rgbindex$PINK) }
            if(aa == "H" & cwcontains(c("%","#","A","C","F","H","I","L","M","V","W","Y","P","p"),consaa)){ return(rgbindex$CYAN) }
            if(aa == "Y" & cwcontains(c("%","#","A","C","F","H","I","L","M","V","W","Y","P","p"),consaa)){ return(rgbindex$CYAN) }
            if(aa == "E" & cwcontains(c("-","E","D","q","Q"),consaa)){ return(rgbindex$MAGENTA) }
            if(aa == "D" & cwcontains(c("-","E","D","n","N"),consaa)){ return(rgbindex$MAGENTA) }
            if(aa == "R" & cwcontains(c("+","K","R","Q"),consaa)){ return(rgbindex$RED) }
            if(aa == "K" & cwcontains(c("+","K","R","Q"),consaa)){ return(rgbindex$RED) }
            return("white")
            return(acol[[aa]])
        }
        '''
    def setupCol(self,overwrite=False): ### Sets up colour lists and dictionaries
        '''Sets up colour lists and dictionaries.'''
        try:### ~ [0] ~ Col list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if overwrite or not self.list['Col']:
                self.list['Col'] = ['rgb(0,0,0)'] + soton[0:]; self.dict['Col']['Col'] = self.list['Col']
                for r in [0,127,255]:
                    for g in [0,127,255]:
                        for b in [0,127,255]:
                             self.list['Col'].append('rgb(%d,%d,%d)' % (r,g,b))
            #!# Add special colour dictionaries #!#
            ### ~ [1] ~ Add ClustalW colours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            acol = self.dict['Col']['AA'] = {}
            '''            
                        for (aa in c("I","L","V","IFLMV","IFLM","IFMV","IFLV","FLMV","ILMV","IL","IV","ILV","ILF")){ acol[aa] = soton$col[19] }
            for (aa in c("A","M")){ acol[aa] = soton$col[19] }
            for (aa in c("AG","GS","AS","AGS")){ acol[aa] = soton$col[21] }
            for (aa in c("K","R","KR")){ acol[aa] = soton$col[3] }
            for (aa in c("D","E","DE")){ acol[aa] = soton$col[15] }
            for (aa in c("S","T","ST")){ acol[aa] = soton$col[14] }
            for (aa in c("C")){ acol[aa] = soton$col[9] }
            for (aa in c("P")){ acol[aa] = soton$col[17] }
            for (aa in c("F","Y","W","FY","FW","WY","FWY")){ acol[aa] = soton$col[18] }
            for (aa in c("G")){ acol[aa] = soton$col[5] }
            for (aa in c("H","HK","HR","HKR")){ acol[aa] = soton$col[2] }
            for (aa in c("HY","FH","FHY")){ acol[aa] = soton$col[1] }
            for (aa in c("Q","N")){ acol[aa] = soton$col[1] }
            for (aa in c("X","x",".")){ acol[aa] = soton$col[6] }
            for (aa in c("0","1","2","3","4","5","6","7","8","9","[","]",",","{","}")){ acol[aa] = "white" }
            for (aa in c("-")){ acol[aa] = "black" }
            for (aa in c("#","^","$","+")){ acol[aa] = soton$col[2] }

            ## ~ Revised CWColour paletter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for (aa in c("I","L","V","IFLMV","IFLM","IFMV","IFLV","FLMV","ILMV","IL","IV","ILV","ILF")){ acol[aa] = cwcol("I",c("#")) }
            for (aa in c("A","M")){ acol[aa] = cwcol("I",c("#")) }
            for (aa in c("AG","GS","AS","AGS")){ acol[aa] = cwcol("A",c("G")) }
            for (aa in c("K","R","KR")){ acol[aa] = cwcol("K",c("+")) }
            for (aa in c("D","E","DE")){ acol[aa] = cwcol("D",c("-")) }
            for (aa in c("S","T","ST")){ acol[aa] = cwcol("S",c("S")) }
            for (aa in c("C")){ acol[aa] = cwcol("C",c("C")) }
            for (aa in c("P")){ acol[aa] = cwcol("P",c("P")) }
            for (aa in c("F","Y","W","FY","FW","WY","FWY")){ acol[aa] = cwcol(substr(aa,1,1),c("#")) }
            for (aa in c("G")){ acol[aa] = cwcol("G",c("G")) }
            for (aa in c("H","HK","HR","HKR")){ acol[aa] = cwcol("H",c("H")) }
            for (aa in c("HY","FH","FHY")){ acol[aa] = cwcol("F",c("#")) }
            for (aa in c("Q","N")){ acol[aa] = cwcol(aa,c(aa)) }

            #ppcol = c("black",soton$col[1:11],soton$col[13:16],soton$col[18:21])
            ppcol = c("black",soton$col[1:11],soton$col[13:21])



            ## ~ ClustalW Alignment colours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # color lookup table - this is optional, if no rgbindex is specified, 8
            # hardcoded colors will be used.
            # A maximum of 16 colors can be specified - any more will be ignored!
            # @rgbindex
            # RED 0.9 0.2 0.1
            # BLUE 0.1 0.5 0.9
            # GREEN 0.1 0.8 0.1
            # CYAN 0.0 0.7 0.7
            # PINK 0.9 0.5 0.5
            # MAGENTA 0.8 0.3 0.8
            # YELLOW 0.8 0.8 0.0
            # ORANGE 0.9 0.6 0.3
            rgbindex = list()
            rgbindex["RED"] = rgb(0.9, 0.2, 0.1)
            rgbindex["BLUE"] = rgb(0.1, 0.5, 0.9)
            rgbindex["GREEN"] = rgb(0.1, 0.8, 0.1)
            rgbindex["CYAN"] = rgb(0.0, 0.7, 0.7)
            rgbindex["PINK"] = rgb(0.9, 0.5, 0.5)
            rgbindex["MAGENTA"] = rgb(0.8, 0.3, 0.8)
            rgbindex["YELLOW"] = rgb(0.8, 0.8, 0.0)
            rgbindex["ORANGE"] = rgb(0.9, 0.6, 0.3)

            # :: @consensus
            cwcons = function(aalist){
                acons = c()
                for(aa in c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")){
                    ax = sum(aalist == aa)
                    acut = 0.85 * length(aalist)
                    if(ax >= acut){ acons = c(acons,aa); }
                }
                nonpol = 0
                for(aa in c("A","C","F","H","I","L","M","P","V","W","Y")){
                    nonpol = nonpol + sum(aalist == aa)
                }
                if(nonpol > 0.6 * length(aalist)){ acons = c(acons,"%"); }
                if(nonpol > 0.8 * length(aalist)){ acons = c(acons,"#"); }
                if(sum(aalist == "E" | aalist == "D") >= (0.5 * length(aalist))) { acons = c(acons,"-"); }
                if(sum(aalist == "K" | aalist == "R") >= (0.6 * length(aalist))) { acons = c(acons,"+"); }
                if(sum(aalist == "G") >= (0.5 * length(aalist))) { acons = c(acons,"g"); }
                if(sum(aalist == "N") >= (0.5 * length(aalist))) { acons = c(acons,"n"); }
                if(sum(aalist == "Q" | aalist == "E") >= (0.5 * length(aalist))) { acons = c(acons,"q"); }
                if(sum(aalist == "P") >= (0.5 * length(aalist))) { acons = c(acons,"p"); }
                if(sum(aalist == "S" | aalist == "T") >= (0.5 * length(aalist))) { acons = c(acons,"s","t"); }
                return(acons)
            }
            # :: % = 60% w:l:v:i:m:a:f:c:y:h:p
            # :: # = 80% w:l:v:i:m:a:f:c:y:h:p
            # :: - = 50% e:d
            # :: + = 60% k:r
            # :: g = 50% g
            # :: n = 50% n
            # :: q = 50% q:e
            # :: p = 50% p
            # :: t = 50% t:s
            # :: s = 50% t:s
            # :: A = 85% a
            # :: C = 85% c
            # :: D = 85% d
            # :: E = 85% e
            # :: F = 85% f
            # :: G = 85% g
            # :: H = 85% h
            # :: I = 85% i
            # :: K = 85% k
            # :: L = 85% l
            # :: M = 85% m
            # :: N = 85% n
            # :: P = 85% p
            # :: Q = 85% q
            # :: R = 85% r
            # :: S = 85% s
            # :: T = 85% t
            # :: V = 85% v
            # :: W = 85% w
            # :: Y = 85% y

            cwcontains = function(wanted,observed){
                for(type in wanted){
                    if(sum(observed == type) > 0){ return(TRUE) }
                }
                return(FALSE)
            }

            # :: @color
            # :: g = ORANGE
            # :: p = YELLOW
            # :: t = GREEN if t:S:T:%:#
            # :: s = GREEN if t:S:T:#
            # :: n = GREEN if n:N:D
            # :: q = GREEN if q:Q:E:+:K:R
            # :: w = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
            # :: l = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
            # :: v = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
            # :: i = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
            # :: m = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
            # :: a = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p:T:S:s:G
            # :: f = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
            # :: c = BLUE if %:#:A:F:H:I:L:M:V:W:Y:S:P:p
            # :: c = PINK if C
            # :: h = CYAN if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
            # :: y = CYAN if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
            # :: e = MAGENTA if -:D:E:q:Q
            # :: d = MAGENTA if -:D:E:n:N
            # :: k = RED if +:K:R:Q
            # :: r = RED if +:K:R:Q
            '''

        except: self.errorLog('Problem during %s setupCol.' % self); return False  # Setup failed
#########################################################################################################################
    def col(self,i=0,type='Col',cycle=True,convert={}):
        try: self.dict['Col'][type][0]
        except: self.setupCol();
        while i >= len(self.dict['Col'][type]): i -= len(self.dict['Col'][type])
        colour = self.dict['Col'][type][i]
        if colour in convert: return convert[colour]
        return colour
#########################################################################################################################
    ### <3> ### SVG Code Methods                                                                                        #
#########################################################################################################################
    def svgHTML(self,svglink,title,svgfile=None,height=1600,width=1600):    ### Returns HTML Code for SVG file embedding.
        '''Returns HTML Code for SVG file embedding.'''
        try:
            if not svgfile: svgfile = svglink 
            svghtm = '<p title="%s">\n' % (title)
            try: (width,height) = rje.matchExp('<svg width="(\d+)" height="(\d+)" version="1.1"',open(svgfile,'r').read())
            except: pass
            svghtm += '<embed src="%s" width="%s" height="%s" type="image/svg+xml"' % (svglink,width,height)
            svghtm += ' pluginspage="http://www.adobe.com/svg/viewer/install/" /></p>'
            return svghtm
        except: self.errorLog(rje_zen.Zen().wisdom()); return '<i>SVG code error!</i>'
#########################################################################################################################
    def svgFile(self,svgtext,filename='',width='100%',height='100%'):  ### Returns SVG code wrapped in header and footer. Saves if filename given.
        '''
        Returns SVG code wrapped in header and footer. Saves if filename given.
        >> svgtext:str = Main body of SVG file
        >> filename:str [''] = Filename to save to. If no filename given, will not save.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            head = ['<?xml version="1.0" standalone="no"?>','<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"',
                    '"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">','',
                    '<svg width="%s" height="%s" version="1.1"' % (width,height),'xmlns="http://www.w3.org/2000/svg">']
            tail = ['','</svg>']
            try: svg = rje.join(head+svgtext+tail,'\n')
            except: svg = rje.join(head+[svgtext]+tail,'\n')
            ### ~ [1] ~ Optionally save file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if filename.lower() not in ['','none']: open(filename,'w').write(svg)
            return svg
        except: self.errorLog('Problem during %s svgFile.' % self)
#########################################################################################################################
    def networkPlot(self,npos,G,nodecol={},font=16,width=1600,height=1200,ntype='ellipse',cutspace=True,xoffset=0,yoffset=0):  ### Plot partial PPI network with given coordinates.
        '''
        Plot partial PPI network with given coordinates.
        >> npos:dict = Dictionary of {node:(x,y)}
        >> G:dict = Dictionary of {node:{node:col}} or {node:[nodes]}
        >> cutspace:bool = whether to cut unneccessary whitespace (crudely)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            svg = ['<!-- SVG Network Plot -->']
            npos = rje.combineDict({},npos)
            font_type = 'Impact'
            ## ~ [0a] ~ Network properties ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            nodes = rje.sortKeys(G)
            spokenum = len(npos)
            xmax = 0.01; ymax = 0.01; nmax = 0
            for xy in npos.values(): xmax = max(xmax,xy[0]); ymax = max(ymax,xy[1])
            for node in npos: nmax = max(nmax,len(node))
            xborder = (nmax * font * 0.5) + 5
            yborder = font + 5
            xscale = (width - 2 * xborder) / xmax; yscale = (height - 2 * yborder) / ymax
            ## ~ [0b] ~ Reduce whitespace regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.deBug(npos)
            if cutspace:    #!# Why does this not work when there are no spokes? #!#
                xlist = []; ylist = []
                for hub in rje.sortKeys(npos):
                    x1 = npos[hub][0]; y1 = npos[hub][1]
                    xlist.append((x1,x1)); ylist.append((y1,y1))
                    if hub not in G: continue
                    i = nodes.index(hub)
                    try: spokes = rje.sortKeys(G[hub]); asdict = True
                    except: spokes = G[hub]; asdict = False
                    for spoke in spokes:
                        if nodes.index(spoke) <= i: continue
                        x2 = npos[spoke][0]; y2 = npos[spoke][1]
                        xlist.append((min(x1,x2),max(x1,x2))); ylist.append((min(y1,y2),max(y1,y2)))
                for zlist in [xlist,ylist]:
                    zlist.sort()
                    #self.deBug(zlist)
                    i = 0
                    while i < len(zlist):
                        while (i+1) < len(zlist) and zlist[i+1][0] <= zlist[i][1]:
                            zlist[i] = (zlist[i][0],max(zlist[i][1],zlist[i+1][1]))
                            zlist.pop(i+1)
                        i += 1
                    #self.deBug(zlist)
                for z in [0,1]:
                    zlist = [xlist,ylist][z]
                    zgap = [1.0/xscale,1.0/yscale][z] * font * 3
                    zshift = {}
                    for i in range(1,len(zlist)):
                        zspace = zlist[i][0]-zlist[i-1][1]
                        if zspace > zgap: zshift[zlist[i-1][1]] = zspace - zgap
                    #self.deBug(zshift)
                    if not zshift: continue
                    for node in npos:
                        znew = zpos = npos[node][z]
                        (xnew,ynew) = npos[node]
                        for zkey in zshift:
                            if zpos > zkey: znew -= zshift[zkey]
                        if z == 0: npos[node] = (znew,ynew)
                        else: npos[node] = (xnew,znew)
                    zshift = sum(zshift.values())
                    if z == 0: xmax -= zshift; xscale = (width - 2 * xborder) / xmax
                    else: ymax -= zshift; yscale = (height - 2 * yborder) / ymax                
            ## ~ [0c] ~ Rescale and add border ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.deBug(npos)
            for node in npos: npos[node] = (npos[node][0]* xscale + xborder + xoffset, npos[node][1] * yscale + yborder + yoffset)
            svg.append('<!-- nmax=%.2f, font=%.1f, width=%d, height=%d -->' % (nmax,font,width,height))
            svg.append('<!-- xmax=%.2f, xscale=%.2f, xborder=%.1f, xoffset=%.1f -->' % (xmax,xscale,xborder,xoffset))
            svg.append('<!-- ymax=%.2f, yscale=%.2f, yborder=%.1f, yoffset=%.1f -->' % (ymax,yscale,yborder,yoffset))
            #self.deBug(rje.join(svg[-3:],'\n'))
            #self.deBug(npos)

            ### ~ [1] ~ Draw connections between nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            svg.append('<!-- Connections between nodes -->')
            for i in range(len(nodes)-1):
                hub = nodes[i]
                x1 = npos[hub][0]; y1 = npos[hub][1]
                try: spokes = rje.sortKeys(G[hub]); asdict = True
                except: spokes = G[hub]; asdict = False
                for spoke in spokes:
                    if nodes.index(spoke) <= i: continue
                    x2 = npos[spoke][0]; y2 = npos[spoke][1]
                    if asdict:
                        if 'style' in G[hub][spoke]: svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" %s/>' % (x1,y1,x2,y2,G[hub][spoke]))
                        else: svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" style="stroke:%s;"/>' % (x1,y1,x2,y2,G[hub][spoke]))
                    else: svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" style="stroke:rgb(0,0,0);stroke-width:2"/>' % (x1,y1,x2,y2))
            svg += ['<!-- END Connections between nodes -->','']

            ### ~ [2] ~ Draw the nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            svg.append('<!-- Draw the nodes -->')
            for node in rje.sortKeys(npos):
                x = npos[node][0]; y = npos[node][1]; ry = font
                if ntype == 'ellipse':
                    rx = max(4,min(7,len(node))) * font * 0.5
                    etxt = '<ellipse cx="%.1f" cy="%.1f" rx="%.1f" ry="%.1f"' % (x,y,rx,ry)
                    if node in nodecol:
                        try: etxt += ' style="fill:%s;' % self.col(int(nodecol[node]),convert={'rgb(0,0,0)':'rgb(255,255,255)'})
                        except: etxt += ' style="fill:%s;' % nodecol[node]
                    else: etxt += ' style="fill:rgb(255,255,255);'
                    etxt += ' stroke:rgb(0,0,0);stroke-width:2"/>'  #!# Could add border colour detail dictionary?
                    svg += ['',etxt]
                if ntype == 'rect':
                    rx = max(4,min(7,len(node))) * font
                    rtxt = '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f"' % (x,y,rx,ry)
                    if node in nodecol: rtxt += ' style="fill:%s;' % nodecol[node]
                    else: rtxt += ' style="fill:rgb(255,255,255);'
                    rtxt += ' stroke:rgb(0,0,0);stroke-width:2"/>'  #!# Could add border
                    svg += ['',rtxt]
                ntxt = '<text x="%.1f" y="%.1f" text-anchor="middle" alignment-baseline="central"' % (x,y)
                if node in nodecol and nodecol[node] in ['rgb(0,0,0)']:
                    ntxt += ' font-family="%s" font-bold="True" font-size="%d" fill="%s" >' % (font_type,font,soton[11])  
                else: ntxt += ' font-family="%s" font-bold="True" font-size="%d" fill="black" >' % (font_type,font)
                ntxt += node
                ntxt += '</text>'
                svg += ['',ntxt]
            svg += ['<!-- END draw nodes -->','']
            return svg
        except: self.errorLog('Problem during %s networkPlot.' % self)
#########################################################################################################################
    def svgTree(self,basefile='',data={},treesplit=0.5,font=12,maxfont=20,width=1600,height=0,save=True,xoffset=0,yoffset=0,internal_labels='boot'):  ### Generate SVG Tree
        '''
        Generate SVG Tree. Based on rje_ppi.r
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            svg = ['<!-- Draw tree -->','']
            ext_font = 'Georgia'; int_font = 'Tahoma'
            nmax = 0
            ### ~ [1] Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup Tree data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not data: data = rje.dataDict(self,'%s.tree.csv' % basefile,['nodenum'],'all')
            ynum = len(data); xmax = 0.0
            for node in rje.sortKeys(data):
                try: int(node)
                except: data.pop(node); continue
                for dkey in data[node]:
                    try:
                        if dkey in ['ypos','xpos','ancy','ancx']: data[node][dkey] = float(data[node][dkey])
                        elif dkey in ['nodenum','anc','boot']: data[node][dkey] = intt(data[node][dkey])
                    except: pass
                xmax = max(xmax,data[node]['xpos'])
                nmax = max(nmax,data[node]['name'])
            xborder = 5; yborder = font + 5
            if not height: height = (font * ynum + 2 * yborder)
            yscale = (height - (2 * yborder)) / float(ynum)        # Add border = font
            xscale = ((treesplit * width) - xborder)  / xmax                 #!# Add border and name space
            xtext = 10; ytext = 0
            maxfont = min(yscale,maxfont)
            xoffset += xborder; yoffset += height
            #!# Add optional Tree title #!#
            ## ~ [1b] Draw tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tree = data
            for i in rje.sortKeys(tree):
                data = tree[i]
                xpos = xoffset + (data['xpos'] * xscale)
                xanc = xoffset + (data['ancx'] * xscale)
                ypos = yoffset - (yborder + (data['ypos'] * yscale) + yscale)
                yanc = yoffset - (yborder + (data['ancy'] * yscale) + yscale)
                
                if xpos != xanc:
                    svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"' % (xpos,ypos,xanc,ypos))
                    if 'col' in data:
                        try: svg.append('style="stroke:%s;stroke-width:2"/>' % (self.col(int(data['col']),convert={'rgb(255,255,255)':'rgb(0,0,0)'})))
                        except: svg.append('style="stroke:%s;stroke-width:2"/>' % (data['col']))
                    else: svg.append('style="stroke:rgb(0,0,0);stroke-width:2"/>')

                svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"' % (xanc,ypos,xanc,yanc))
                if 'col' in data:
                    try: svg.append('style="stroke:%s;stroke-width:2"/>' % (self.col(int(data['col']),convert={'rgb(255,255,255)':'rgb(0,0,0)'})))
                    except: svg.append('style="stroke:%s;stroke-width:2"/>' % (data['col']))
                else: svg.append('style="stroke:rgb(0,0,0);stroke-width:2"/>')
            ## ~ [1c] Add text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if int(i) <= (ynum+1) / 2.0: 
                    ntxt = '<text x="%.1f" y="%.1f" text-anchor="left" alignment-baseline="central"' % (xpos+xtext,ypos+ytext)
                    if 'col' in data:
                        try: ntxt += ' font-family="%s" font-bold="True" font-size="%d" fill="%s" >' % (ext_font,max(maxfont,font),self.col(int(data['col']),convert={'rgb(255,255,255)':'rgb(0,0,0)'}))
                        except: ntxt += ' font-family="%s" font-bold="True" font-size="%d" fill="%s" >' % (ext_font,max(maxfont,font),data['col'])
                        try:    # Update ancestral colours #
                            anci = int(data['anc'])
                            if 'col' in tree[anci] and tree[anci]['col'] and tree[anci]['col'] != data['col']: tree[anci]['col'] = 'black'
                            else: tree[anci]['col'] = data['col']
                        except: pass
                    else: ntxt += ' font-family="%s" font-bold="True" font-size="%d" fill="black" >' % (ext_font,max(maxfont,font))
                    ntxt += data['name']
                    ntxt += '</text>'
                    svg.append(ntxt)
                else:
                    ntxt = '<text x="%.1f" y="%.1f" text-anchor="right" alignment-baseline="central"' % (xpos+xtext,ypos+ytext)
                    if 'col' in data:
                        try: ntxt += ' font-family="%s" font-bold="True" font-style="italic" font-size="%d" fill="%s" >' % (int_font,max(maxfont,font)*0.7,self.col(int(data['col'])))
                        except: ntxt += ' font-family="%s" font-bold="True" font-style="italic" font-size="%d" fill="%s" >' % (int_font,max(maxfont,font)*0.7,data['col'])  
                    else: ntxt += ' font-family="%s" font-bold="True" font-style="italic" font-size="%d" fill="%s" >' % (int_font,max(maxfont,font)*0.7,self.col(5))
                    try: svg.append(ntxt + data[internal_labels] + '</text>')
                    except: self.errorLog('Internal "%s" label error!' % internal_labels)
            ## ~ [1d] Add scale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            xbars = [0.2,0.5,1.0,2.0,5.0,10.0,25.0,50.0,100.0]
            xbar = 0.1
            while xbar * xscale < 2 * font and xbars: xbar = xbars.pop(0)
            svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"' % (xoffset,yoffset - yborder*2,xoffset+xbar*xscale,yoffset - yborder*2))
            svg.append('style="stroke:%s;stroke-width:2"/>' % (soton[6]))
            svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"' % (xoffset,yoffset - yborder*2,xoffset,yoffset - yborder))
            svg.append('style="stroke:%s;stroke-width:2"/>' % (soton[6]))
            svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"' % (xoffset+xbar*xscale,yoffset - yborder*2,xoffset+xbar*xscale,yoffset - yborder))
            svg.append('style="stroke:%s;stroke-width:2"/>' % (soton[6]))
            ntxt = '<text x="%.1f" y="%.1f" text-anchor="middle" alignment-baseline="central"' % (xoffset,yoffset - yborder/2)
            ntxt += ' font-family="%s" font-bold="True" font-size="%d" fill="%s" >0</text>' % (int_font,font,soton[6])
            svg.append(ntxt)
            ntxt = '<text x="%.1f" y="%.1f" text-anchor="middle" alignment-baseline="central"' % (xoffset+xbar*xscale,yoffset - yborder/2)
            ntxt += ' font-family="%s" font-bold="True" font-size="%d" fill="%s" >%s</text>' % (int_font,font,soton[6],xbar)
            svg.append(ntxt)
            svg += ['<!-- END draw tree -->','']
            if save: self.svgFile(rje.join(svg,'\n'),filename='%s.svg' % basefile,width=width,height=height)
            return svg
        except: self.errorLog('Error in %s.svgTree' % self)
#########################################################################################################################
    def svgAlignment(self): ###
        '''.'''
        return
#########################################################################################################################
### End of SECTION II: SVG Class                                                                                        #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
'''
## ~ Generate sequence logo from profile data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
seqLogo = function(infile,fex=1){
	## ~ Read Data and Setup ~ ##
	motif = read.table(infile,sep="\t",header=TRUE)
	aas = colnames(motif)[-1]
	x = 1.0 / length(motif[,1])
	## ~ Create empty plot ~ ##
	old_mar = par()$mar
	par(mar=c(0,2,2,0))
	#!# png(filename = outfile, width=50*length(motif[,1]), height=500, units = "px")	#, pointsize=25)
	plot(0:1,0:1,type="n",axes=FALSE,ann=FALSE,mar=c(0.1,0.5,0.5,0.1))
	axis(side=2,lwd=2,cex=fex*2)
	## ~ Plot Profile ~ ##
	ix = 0.0
	for (i in 1:length(motif[,1])){
		## Normalise data ##
		isum = sum(motif[i,-1])
		if (isum > 0){ 
            	for (a in aas){ motif[i,a] = motif[i,a] / isum }
        	}
		## Plot ##
		iy = 0.0
		for (a in order(rank(motif[i,]))){
			aa = colnames(motif)[a]
			if (aa == "Pos"){
				text(ix+(x/2),1.05,motif[i,1],adj=c(0.5,0.5),xpd=TRUE,font=2,family="mono",cex=fex*min(1,4/nchar(motif[i,1])))
			}else{
				f = motif[i,aa]
				rect(ix,iy,ix+x,iy+f,col=acol[[aa]],lwd=2)
				if (f > 0.05){ text(ix+(x/2),iy+(f/2),aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex) }
				iy = iy + f
			}
		}
		ix = ix + x
	}
	#!# dev.off()
	par(mar=old_mar)
	length(motif[,1])
}

## ~ Draw a wrapped multiple sequence alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
wrapAln = function(infile,wrap=80,fex=1,start=1,end=0,usecwcol=FALSE){
	## Load Data and Setup ##
	aln = read.table(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	seqnum = length(aln[1,])
	seqlen = length(aln[,1])
	if (wrap<1){ wrap = seqlen - start + 1 }
	if (end<start){ end = seqlen }
	blocks = as.integer((end-start)/wrap) + 1
	height = (blocks * (seqnum + 2)) - 1
	## Create empty plot ##
	old_mar = par()$mar
	par(mar=c(0,5,1.1,0))
	#!# png(filename = outfile, width=20*wrap, height=height*30, units = "px")
	par(mar=c(0,5,1.1,0))
	plot(c(0,wrap),c(0,height),type="n",axes=FALSE,ann=FALSE)
	## Plot alignment ##
	for (b in 1:blocks){
		h = height - ((b-1)*(seqnum+2)) - 1
		for (a in 1:seqnum){
			name = colnames(aln)[a]
			text(-0.5,h-a+0.5,colnames(aln)[a],adj=c(1,0.5),xpd=TRUE,font=1,family="sans",cex=fex*min(1,12/seqnum))
			for (i in 1:wrap){
				p = (b-1)*wrap + i + start - 1
				if (p > seqlen){ break }
				if (aln[p,a] == 0 && substr(name,nchar(name)-3,nchar(name)) == ".ppi"){ aa = "-" }
				else { aa = as.character(aln[p,a]) }
				#if(usecwcol){ rect(i-1,h-a,i,h-a+1,col=cwcol(aa,cwcons(aln[p,2:seqnum])),border=NA) }
				#if(usecwcol){ rect(i-1,h-a,i,h-a+1,col=cwcol(aa,c(aa)),border=NA) }
				#else{ rect(i-1,h-a,i,h-a+1,col=acol[[aa]],border=NA) }
				rect(i-1,h-a,i,h-a+1,col=acol[[aa]],border=NA)
				if (aa == "-"){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",col="white",cex=fex)
				}
				if (aa != "-" && nchar(aa) <= 1){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex)
				}
				if (aa != "-" && nchar(aa) > 1){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",srt=90,cex=fex*2/nchar(aa))
				}
				if (i == 1 | i == wrap | as.integer(p/10) == p/10 | p == seqlen){ text(i-0.5,h+0.7,p,adj=c(0.5,1),xpd=TRUE,family="sans",cex=fex*0.8) }
			}
		}
	}
	#!# dev.off()
	par(mar=old_mar)
	height
}
#x#alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=80)

##### ~~~ New wrapped Alignment function ~~~~~~~~~~~~~~~~~~~~~~~~~~ #####
wrapAlnPNG = function(basefile,wrap=80,mode="png",usecwcol=FALSE){
	if (mode=="test"){
		wrap=120
		basefile = "SLiMJIM/html/visualisations/PPP1CA.TAI12_HUMAN__Q9H175"
	}
	infile = paste(basefile,".aln.tdt",sep="")
	outfile = paste(basefile,".png",sep="")
	if (mode=="test"){ outfile = "test.png"; mode="png" }
	### ~ [1] ~ Load Data and Setup ~~~~~~~~~~~~~~~~~~~~~~~~~ ###
	aln = read.table(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	## ~ [1a] ~ Alignment stats ~~~~~~~~~~~~~~~~~~ ##
	seqnum = length(aln[1,])				# Number of sequences
	seqlen = length(aln[,1])				# Length of sequences
	namelen = 0							# Length of longest name
	for (a in 1:seqnum){ namelen = max(strwidth(colnames(aln)[a],units="inches"),namelen) }
	if (wrap == 0){ wrap = seqlen }			# No wrap = full length
	blocks = as.integer((seqlen-1)/wrap) + 1		# No of wrapped alignment blocks
	height = (blocks * (seqnum + 2)) - 1		# Height (alignment lines) of wrapped alignment
	## ~ [1b] ~ Calculate plot characteristics ~~~ ##
	fex = 0.8 #0.16 / par()$csi		# Should be 0.8 normally?
	xin = wrap / 8				# Width of plot (inches)
	yin = height / 5				# Height of plot (inches)
	namelen = fex * (namelen + 0.0625)	# Length needed for sequence names (inches)
	if (mode=="size"){ return(c(xin,yin)) }
	## Create empty plot ##
	old_par = par()
	#par(omi=c(0,0,0,0))
	par(mar=c(0.1,5,0.1,0.1))
	if (mode=="png"){ png(filename=outfile, width=xin, height=yin, units="in", res=150) }
	if (mode=="test"){ 
		png(filename=outfile, width=xin, height=yin, units="in", res=150) 
	}
	par(mar=c(0.1,5,0.1,0.1))
	plot(c(0,wrap),c(0,height),type="n",axes=FALSE,ann=FALSE,asp=yin/xin)	#pin=c(xin,yin))
	## Plot alignment ##
	for (b in 1:blocks){
		h = height - ((b-1)*(seqnum+2)) - 1
		for (a in 1:seqnum){
			name = colnames(aln)[a]
			text(-0.5,h-a+0.5,name,adj=c(1,0.5),xpd=TRUE,font=1,family="sans",cex=fex)	#min(1,12/seqnum))
			for (i in 1:wrap){
				p = (b-1)*wrap + i
				if (p > seqlen){ break }
				if (aln[p,a] == 0 && substr(name,nchar(name)-3,nchar(name)) == ".ppi"){ aa = "-" }
				else { aa = as.character(aln[p,a]) }
				rect(i-1,h-a,i,h-a+1,col=acol[[aa]],border=NA)
				if (aa == "-"){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",col="white",cex=fex)
				}
				if (aa != "-" && nchar(aa) <= 1){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex)
				}
				if (aa != "-" && nchar(aa) > 1){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",srt=90,cex=2*fex/nchar(aa))
				}
				if (i == 1 | i == wrap | as.integer(p/10) == p/10 | p == seqlen){ text(i-0.5,h+0.5,p,adj=c(0.5,0.5),xpd=TRUE,family="sans",cex=fex) }
			}
		}
	}
	if (mode=="png"){ 
		dev.off() 
	}
	par(mar=old_par$mar,omi=old_par$omi)
	height
}

########## ~~~ MOTIF MULTIPANEL PLOT ~~ ##############
multiPanelMotif = function(basefile){
	## ~ Setup multipanel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	#alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=80)
	pngfile=paste(basefile,".png",sep="")
    #CairoPNG(filename=pngfile, width=2400, height=1600, pointsize=12)
	png(filename=paste(basefile,".png",sep=""), width=2400, height=1600, units = "px", pointsize=12)
	#panels = c(1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,3,3,3,4,4,4,5,5,5,3,3,3,4,4,4,5,5,5)
	panels = c(1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,2,2,2,3,4,4,4,5,5,5,5,5,3,4,4,4,5,5,5,5,5,3,4,4,4,5,5,5,5,5)
	layout(matrix(panels,byrow=TRUE,nrow=5))
	#layout.show(11)
	fex = 2.67

	## ~ Plot alignment surrounding motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=0,fex=fex)

	## ~ Plot sequence logo of motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	mlen = seqLogo(paste(basefile,".profile.tdt",sep=""),fex=fex)

	## ~ Plot tree for motif-containing spokes (based on GABLAM) ~~~~~~~~~ ##
	makeTree(paste(basefile,".tree.csv",sep=""),fex=fex)

	## ~ Make Heatmap (based on GABLAM) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	heatMap(paste(basefile,".heatmap.tdt",sep=""),fex=fex)

	## ~ Plot PPI network for motif-containing spokes ~~~~~~~~~~~~~~~~~~~~ ##
	networkRing(paste(basefile,".ppi.tdt",sep=""),fex=fex)
	dev.off()
}
if (slimjim == "motif"){ multiPanelMotif(basefile) }

######## ~~~ RELATIVE CONSERVATION AND DISORDER PLOT ~~~~ #########################
relConsPlot = function(infile,start=1,end=0,fex=1){
	## Load Data and Setup ##
	rel = read.table(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE,as.is=TRUE)
	seqlen = length(rel[,1])
	if (end<1){ end = seqlen }
	wrap = end - start
	## Create empty plot ##
	old_mar = par()$mar
	par(mar=c(0,5,1.1,0))
	#!# png(filename = outfile, width=20*wrap, height=height*30, units = "px")
	par(mar=c(0,5,1.1,0))
	plot(c(0,wrap),c(-1.5,1.5),type="n",axes=FALSE,ann=FALSE)
	axis(side=2,lwd=2,cex=fex)
	## Plot alignment ##
	for (i in 0:wrap){
		r = start + i
		g = 0.1	# Size of gap each side of bar
		dis = rel[r,"IUPred"]
		ord = rel[r,"Disorder"]
		pos = rel[r,"Pos"]
		aa = rel[r,2]
		cons = rel[r,"RelVNE"]
		# Plot relcons, then disorder, then sequence on top. #
		if (cons<0){
			rect(i+g,0,i+1-g,cons,lwd=2,border=soton$col[5],col=soton$col[9])
		}else{
			rect(i+g,0,i+1-g,cons,lwd=2,border=soton$col[5],col=soton$col[5])
		}
		if (ord == "Order"){
			rect(i+g,0,i+1-g,dis,lwd=2,border=soton$col[1],col=soton$col[2])
		}else{
			rect(i+g,0,i+1-g,dis,lwd=2,border=soton$col[1])
		}
		if (cons>0 & cons<dis & ord=="Order"){ rect(i+g,0,i+1-g,cons,lwd=2,border=soton$col[5]) }
		# Plot sequence and position #
		rect(i,-0.01,i+1,-0.11,col="white",border="white")
		rect(i,-0.015,i+1,-0.105,col=acol[[aa]],border=NA)
		text(i+0.5,-0.06,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex)
		if (i == 0 | i == wrap | as.integer(r/10) == r/10 | r == seqlen){ text(i+0.5,-0.16,r,adj=c(0.5,0.5),xpd=TRUE,family="sans",cex=fex) }
	}
	par(mar=old_mar)
}

######## ~~~ RELATIVE CONSERVATION, DISORDER AND ALIGNMENT PLOT ~~~~ #########################
relConsAln = function(basefile,start=1,end=0,fex=1,lex=2){
	## Load Data and Setup ##
	rel = read.table(paste(basefile,".rel.tdt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE,as.is=TRUE)
	aln = read.table(paste(basefile,".aln.tdt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE,as.is=TRUE)
	seqlen = length(rel[,1])
	if (end<1){ end = seqlen }
	wrap = end - start
	## Create empty plot ##
	old_mar = par()$mar
	par(mar=c(0,5,1.1,0),ljoin="mitre")
	#!# png(filename = outfile, width=20*wrap, height=height*30, units = "px")
	par(mar=c(0,5,1.1,0))
	plot(c(0,wrap),c(-1.5,1.5),type="n",axes=FALSE,ann=FALSE)
	g = 0.1	# Size of gap each side of bar
	axis(side=2,lwd=2,cex=fex,pos=-(2*g))
        rect(-(4*g)-1,-1.6,-(4*g),1.6,col="white",xpd=TRUE,border=NA)
 	for (y in c(1.5,1.0,0.5,0.0,-0.5,-1.0,-1.5)){
 	        text(-(5*g),y,y,adj=c(1,0.5),xpd=TRUE,family="sans",cex=0.8*fex)
        }
	## Plot alignment ##
	for (i in 0:wrap){
		r = start + i
		if(r>seqlen){ break }
		dis = rel[r,"IUPred"]
		ord = rel[r,"Disorder"]
		pos = rel[r,"Pos"]
		aa = rel[r,2]
		cons = rel[r,"RelVNE"]
		# Plot relcons, then disorder, then sequence on top. #
		if (cons<0){
			rect(i+g,0,i+1-g,cons,lwd=lex,border=soton$col[5],col=soton$col[9])
		}else{
			rect(i+g,0,i+1-g,cons,lwd=lex,border=soton$col[5],col=soton$col[5])
		}
		if (ord == "Order"){
			rect(i+g,0,i+1-g,dis,lwd=lex,border=soton$col[1],col=soton$col[2])
		}else{
			rect(i+g,0,i+1-g,dis,lwd=lex,border=soton$col[1])
		}
		if (cons>0 & cons<dis & ord=="Order"){ rect(i+g,0,i+1-g,cons,lwd=lex,border=soton$col[5]) }
		# Plot sequence and position #
		rect(i,-0.01,i+1,-0.11,col="white",border="white")
		rect(i,-0.015,i+1,-0.105,col=acol[[aa]],border=NA)
		text(i+0.5,-0.06,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex)
		if (i == 0 | i == wrap | as.integer(r/10) == r/10 | r == seqlen){ text(i+0.5,-0.16,r,adj=c(0.5,0.5),xpd=TRUE,family="sans",cex=fex) }
		# Add Motif #
		if (aln[r,1] != "-" & as.character(aln[r,1]) != "0"){
            	h = min(1.5,max(cons,dis) + 0.05)
			rect(i,h,i+1,h+0.09,col=acol[[aln[r,1]]],border=NA,xpd=TRUE)
		   	text(i+0.5,h+0.045,aln[r,1],adj=c(0.5,0.5),xpd=TRUE,font=2,family="mono",cex=fex)
            }

	}
	## Add key ##
	k = 0.02
	rect(k,-1.0-k,11+k,-1.125+k,lwd=2,border=NA,col="white")
	rect(k,-1.0-k,1+k,-1.125+k,lwd=2,border=soton$col[5],col=soton$col[5])
	text(1.5+k,-1.0625,"Conserved",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	rect(k,-1.125-k,11+k,-1.25+k,lwd=2,border=NA,col="white")
	rect(k,-1.125-k,1+k,-1.25+k,lwd=2,border=soton$col[5],col=soton$col[9])
	text(1.5+k,-1.1875,"Unconserved",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	rect(k,-1.25-k,11+k,-1.375+k,lwd=2,border=NA,col="white")
	rect(k,-1.25-k,1+k,-1.375+k,lwd=2,border=soton$col[1],col=soton$col[2])
	text(1.5+k,-1.3125,"Ordered",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	rect(k,-1.375-k,11+k,-1.5+k,lwd=2,border=NA,col="white")
	rect(k,-1.375-k,1+k,-1.5+k,lwd=2,border=soton$col[1])
	text(1.5+k,-1.4375,"Disordered",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	par(mar=old_mar)
}

########## ~~~ SPOKE ALIGN PLOT ~~ ##############

testSpokeAln = function(basefile){
	## ~ Single panel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	aln = read.table(paste(basefile,".aln.tdt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE)
	seqnum = length(aln[1,])
	seqlen = length(aln[,1])
	wrap = 80
	if (wrap<1){ wrap = seqlen }
	blocks = as.integer((seqlen-1)/wrap) + 1
	alnh = (blocks * (seqnum + 2)) - 1
	png(filename=paste(basefile,".png",sep=""), width=20+20*wrap, height=20*alnh, units = "px", pointsize=12)
	alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=wrap,fex=1.5)
	dev.off()
	png(filename=paste(basefile,".rel.png",sep=""), width=20+20*seqlen, height=600, units = "px", pointsize=12)
	relConsAln(basefile,fex=0.8,lex=2)
	dev.off()


	## ~ Setup multipanel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	wrap = 80
	if (wrap<1){ wrap = seqlen }
	blocks = as.integer((seqlen-1)/wrap) + 1
	alnh = (blocks * (seqnum + 2)) - 1
	png(filename=paste(basefile,".2.png",sep=""), width=20+20*wrap, height=20*alnh+400*blocks, units = "px", pointsize=12)
	panels = c(rep(1,blocks),2:(blocks+1))
	layout(matrix(panels,byrow=TRUE,ncol=1))
	#par(mfcol=c(blocks*2,1))
	wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=wrap,fex=1.5)
	for (b in 1:blocks){
		relConsAln(basefile,start=1+(b-1)*wrap,end=b*wrap)	#!# Add splitting and wrapping #!#
	}
	dev.off()

	## ~ Setup multipanel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	wrap = 0
	if (wrap<1){ wrap = seqlen }
	blocks = as.integer((seqlen-1)/wrap) + 1
	alnh = (blocks * (seqnum + 2)) - 1
	png(filename=paste(basefile,".3.png",sep=""), width=20+20*wrap, height=20*alnh+400*blocks, units = "px", pointsize=12)
	panels = c(1:2)
	layout(matrix(panels,byrow=TRUE,ncol=1))
	#par(mfcol=c(blocks*2,1))
	wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=wrap,fex=1.5)
	for (b in 1:blocks){
		relConsAln(basefile,start=1+(b-1)*wrap,end=b*wrap)	#!# Add splitting and wrapping #!#
	}
	dev.off()
}

spokeAln = function(basefile){
	## ~ Single panel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	aln = read.table(paste(basefile,".aln.tdt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE)
	seqnum = length(aln[1,])
	seqlen = length(aln[,1])
	start = 1
	end = 110
	wrap = 110
	while (start < seqlen){
		pngfile = paste(basefile,".",preZero(start),"-",preZero(end),".png",sep="")
		alnh = seqnum + 2
        #CairoPNG(filename=pngfile, width=2400, height=min(1200,80*alnh), pointsize=12)
        png(filename=pngfile, width=2400, height=min(1200,80*alnh), units = "px", pointsize=12)
		panels = c(rep(1,as.integer((seqnum-1)/25)+1),rep(2,2))
		layout(matrix(panels,byrow=TRUE,ncol=1))
		wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=wrap,fex=3,start=start,end=end,usecwcol=TRUE)
		relConsAln(basefile,fex=3,lex=2,start=start,end=end)
		dev.off()
		if (end > seqlen){ break }
		start = start + 100
		end = end + 100
	}
}
if (slimjim == "spokealn"){ spokeAln(basefile) }

'''
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: print('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:
        print('\n\n *** No standalone functionality! *** \n\n')

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print('Cataclysmic run error:', sys.exc_info()[0])
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
