#!/usr/local/bin/python

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
# Author contact: <redwards@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       RJE_XGMML
Description:  RJE XGMLL Module 
Version:      1.0
Last Edit:    26/06/14
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is currently designed to store data for, and then output, an XGMML file for uploading into Cytoscape etc.
    Future versions may incoporate the ability to read and manipulate existing XGMML files.

Commandline:
    At present, all commands are handling by the class populating the XGMML object.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, math, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Basic functional version for use with other modules. Disabled attributes by default for Cytoscape.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Read XGMML back into object
    # [ ] : Read and convert other formats
    # [ ] : Read distance matrices, perform MDS and output?
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyyear) = ('RJE_XGMML', '1.0', 'June 2014', '2007')
    description = 'RJE XGMLL Module'
    author = 'Dr Richard J. Edwards.'
    comments = [rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('{0} v{1}'.format(info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: XGMML Class                                                                                             #
#########################################################################################################################
class XGMML(rje.RJE_Object):     
    '''
    XGMML Class. Author: Rich Edwards (2007).

    Info:str
    - Name = This is the ID used for the graph [RJE_XGMML]
    - Description = Description of network
    - Type = Type of data in network
    
    Opt:boolean
    - XGMMLAtt = Whether to output full XGMML attributes, e.g. colour/size etc. [False]

    Stat:numeric

    List:list

    Dict:dictionary
    - Edge = Dictionary of edges between nodes {Type:{(source,target):Attributes}}
    - EdgeAtt = Dictionary of edge attributes {Att:Type}
    - Node = Dictionary of Nodes to be output {Node:Attributes}
    - NodeAtt = Dictionary of node attributes {Att:Type}
    - NodePos = Dictionary of node positions {Node:[x,y]}

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Description','Type']
        self.optlist = ['XGMMLAtt']
        self.statlist = []
        self.listlist = []
        self.dictlist = ['Edge','EdgeAtt','Node','NodeAtt','NodePos']
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Name':'RJE_XGMML','Description':'RJE XGMML output for Cytoscape','Type':'Mixed network data'})
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ###
                self._cmdReadList(cmd,'opt',['XGMMLAtt'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Output Methods                                                                                          #
#########################################################################################################################
    def saveXGMML(self,filename=None,format='Cytoscape'):       ### Saves object data to file in XGMML format
        '''
        Saves object data to file in XGMML format.
        >> filename:str [None] = Output file. Will use name.xgmml if None.
        >> format:str [Cytoscape] = Target for output file
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not filename or filename.lower() == 'none': filename = '%s.xgmml' % self.info['Name']
            self.log.printLog('#XGMML','Output of XGMML file %s for %s...' % (filename,format),log=False,newline=False)
            rje.backup(self,filename)
            date = rje.dateTime()
            OUT = open(filename,'w')
            
            ### ~ [2] Output headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            OUT.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
            OUT.write('<graph label="%s" id="%s" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns="http://www.cs.rpi.edu/XGMML" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n' % (self.info['Name'],self.info['Name']))
            OUT.write('    <att name="documentVersion" value="1.0"/>\n')
            ## ~ [2a] Cytoscape format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            OUT.write('    <att name="networkMetadata">\n')
            OUT.write('        <rdf:RDF>\n')
            OUT.write('            <rdf:Description rdf:about="http://www.cytoscape.org/">\n')
            OUT.write('                <dc:source>RJE_XGMML</dc:source>\n')
            OUT.write('                <dc:format>Cytoscape-XGMML</dc:format>\n')
            OUT.write('                <dc:description>%s</dc:description>\n' % self.info['Description'])  
            OUT.write('                <dc:date>%s</dc:date>\n' % date)
            OUT.write('                <dc:type>%s</dc:type>\n' % self.info['Type'])
            OUT.write('                <dc:identifier>N/A</dc:identifier>\n')
            OUT.write('                <dc:title>%s</dc:title>\n' % self.info['Name'])
            OUT.write('            </rdf:Description>\n')
            OUT.write('        </rdf:RDF>\n')
            OUT.write('    </att>\n\n')

            ### ~ [3] Output Nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            size = 35.0
            nodelist = rje.sortKeys(self.dict['Node'])
            (n,x,y) = (int(math.sqrt(len(nodelist))),0,0)
            for node in nodelist:
                try:
                    ## ~ [3a] Basic node attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    OUT.write('    <node label="%s" id="%d">\n' % (node,nodelist.index(node)))
                except: self.errorLog('!'); continue
                try:
                    for att in rje.sortKeys(self.dict['Node'][node]):
                        if att not in self.dict['NodeAtt']: continue
                        type = self.dict['NodeAtt'][att]
                        value = rje.replace('%s' % self.dict['Node'][node][att],'&','and')
                        OUT.write('        <att type="%s" name="%s" label="%s" value="%s"/>\n' % (type,att,att,value))
                    ### ~ [3b] Graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    #!# Add control for these at some point! #!#
                    if node in self.dict['NodePos']:
                        (nx,ny) = self.dict['NodePos'][node]
                        try: nx * size
                        except: self.errorLog('%s nodepos X = %s' % (node,nx)); nx = x
                        try: ny * size
                        except: self.errorLog('%s nodepos Y = %s' % (node,ny)); ny = y
                    else: [nx,ny] = [x,y]
                    if self.getBool('XGMMLAtt'):
                        OUT.write('        <graphics w="%.1f" h="%.1f" width="1" type="ellipse" outline="#000000" fill="#ff9999" y="%.1f" x="%.1f">\n' % (size,size,ny*2*size,nx*2*size))
                        OUT.write('            <att name="cytoscapeNodeGraphicsAttributes">\n')
                        OUT.write('                <att name="nodeTransparency" value="1.0"/>\n')
                        #OUT.write('                <att name="nodeLabelFont" value="Default-0-12"/>\n')
                        OUT.write('                <att name="borderLineType" value="solid"/>\n')
                        OUT.write('            </att>\n')
                    else:
                        OUT.write('        <graphics y="%.1f" x="%.1f">\n' % (ny*2*size,nx*2*size))
                    OUT.write('        </graphics>\n')
                    x += 1
                    if x > n: (x,y) = (0,y+1)
                    ### ~ [3c] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                except: self.errorLog('!')
                OUT.write('    </node>\n')
            
            ### ~ [4] Output Edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for etype in rje.sortKeys(self.dict['Edge']):
                for edge in rje.sortKeys(self.dict['Edge'][etype]):
                    try:
                        ## ~ [3a] Basic edge attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        id = '%s (%s) %s' % (edge[0],etype,edge[1])
                        OUT.write('    <edge label="%s" id="%s" target="%d" source="%d">\n' % (id,id,nodelist.index(edge[1]),nodelist.index(edge[0])))
                    except: self.errorLog('!'); continue
                    try:
                        OUT.write('        <att type="string" name="canonicalName" label="canonicalName" value="%s"/>\n' % id)
                        OUT.write('        <att type="string" name="TYPE" label="TYPE" value="%s"/>\n' % etype)
                        if 'interaction' not in self.dict['EdgeAtt']: OUT.write('        <att type="string" name="interaction" label="interaction" value="%s"/>\n' % etype)
                        OUT.write('        <att type="string" name="EDGE_TYPE" label="EDGE_TYPE" value="DefaultEdge"/>\n')
                        for att in self.dict['Edge'][etype][edge]:
                            if att.lower() == 'type': continue
                            if att not in self.dict['EdgeAtt']: continue
                            type = self.dict['EdgeAtt'][att]
                            value = rje.replace('%s' % self.dict['Edge'][etype][edge][att],'&','and')
                            OUT.write('        <att type="%s" name="%s" label="%s" value="%s"/>\n' % (type,att,att,value))
                        ### ~ [3b] Graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        #!# Update these at some point! #!#
                        #OUT.write('        <graphics width="1" fill="#0000ff">\n')
                        #OUT.write('            <att name="cytoscapeEdgeGraphicsAttributes">\n')
                        #OUT.write('                <att name="sourceArrow" value="0"/>\n')
                        #OUT.write('                <att name="targetArrow" value="0"/>\n')
                        #OUT.write('                <att name="edgeLabelFont" value="Default-0-10"/>\n')
                        #OUT.write('                <att name="edgeLineType" value="SOLID"/>\n')
                        #OUT.write('                <att name="sourceArrowColor" value="#000000"/>\n')
                        #OUT.write('                <att name="targetArrowColor" value="#000000"/>\n')
                        #OUT.write('                <att name="curved" value="STRAIGHT_LINES"/>\n')
                        #OUT.write('            </att>\n')
                        #OUT.write('        </graphics>\n')
                        ### ~ [3c] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    except: self.errorLog('!')
                    OUT.write('    </edge>\n')
                            
            ### ~ [5] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            OUT.write('</graph>\n')
            OUT.close()
            self.log.printLog('\r#XGMML','Output of XGMML file %s for %s complete.' % (filename,format))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: XGMML Class                                                                                      #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return

    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: print('\n\n *** No standalone functionality! *** \n\n')
        
    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
