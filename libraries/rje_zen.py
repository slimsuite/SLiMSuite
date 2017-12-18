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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 31 Shanagarry, Milltown Road, Milltown, Dublin 6, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_zen
Description:  Random Zen Wisdom Generator
Version:      1.3.2
Last Edit:    13/07/17
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    Generates random (probably nonsensical) Zen wisdoms. Just for fun.

Commandline:
    wisdoms=X   : Number of Zen Wisdoms to return [10]
    zensleep=X  : Time in seconds to sleep between wisdoms [0]

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, random, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Full working version with four zen types and reasonable vocabulary.
    # 1.1 - Added a few more words here and there.
    # 1.2 - Added a webserver mode to return text directly.
    # 1.3.0 - Modified output to work with new REST service calls.
    # 1.3.1 - Added some more words.
    # 1.3.2 - Added some more words.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add reading of vocabulary from files?
    # [ ] : Add "together, they fight crime!" type of messages.
    # [ ] : Document constructions better.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('RJE_ZEN', '1.3.2', 'July 2017', '2007')
    description = 'Random Zen Wisdom Generator'
    author = 'Dr Richard J. Edwards.'
    comments = ['WARNING: These wisdoms are computer-generated garbage.', 'Heed them at your own peril.']
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
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
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Zen Class                                                                                               #
#########################################################################################################################
class Zen(rje.RJE_Object):     
    '''
    Zen Wisdom Class. Author: Rich Edwards (2005).

    Info:str
    
    Opt:boolean

    Stat:numeric
    - Wisdoms = Number of Zen Wisdoms to return [10]
    - ZenSleep = Time in seconds to sleep between wisdoms [0]

    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = []
        self.optlist = []
        self.statlist = ['Wisdoms','ZenSleep']
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setStat({'Wisdoms':10,'ZenSleep':0})
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)
                self._cmdReadList(cmd,'int',['Wisdoms','ZenSleep'])
            except:self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text.
        Individual outputs can be identified/parsed:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ...

        Outputs available:
            wisdoms = the list of randomly generated nonsense.

        &rest=OUTFMT can then be used to retrieve individual parts of the output in future.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.info key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    ### <2> ### Main Run Methods                                                                                        #
#########################################################################################################################
    def run(self):    ### Main run method
        '''Main Run method. Calls wisdom method X times.'''
        try:### ~ Call and print wisdoms ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            wisdoms = []
            for w in range(self.stat['Wisdoms']):
                zenwisdom = self.wisdom()
                wisdoms.append('%s.' % zenwisdom)
                if self.getStrLC('Rest') or not self.server(): self.log.printLog('#ZEN',zenwisdom)
                time.sleep(self.getInt('ZenSleep'))
            self.dict['Output']['wisdoms'] = string.join(wisdoms,'\n')
            if self.server(): return string.join(wisdoms,'\n')
            if self.stat['Interactive'] >= 0 and self.opt['Win32']: rje.choice('\n<ENTER> to Quit')
        except: self.log.errorLog('Bad vibes from Zen.run()')
#########################################################################################################################
    def wisdom(self,wtype=None):      ### Generates and returns a random Zen wisdom
        '''Generates and returns a random Zen wisdom'''
        try:### ~ Call the appropriate Zen Method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if type in ['A','B','C','D']: mytype = wtype
            else: mytype = rje.randomList(['A'] * 6 + ['B'] * 7 + ['C'] * 3 + ['D'] * 1)[0]
            if mytype == 'A': zen = self._zenA()    # Type A = 'The WISE MAN X BUT THE WWW XXX Y'  
            if mytype == 'B': zen = self._zenB()    # Type B = 'It is better to X when Y'
            if mytype == 'C': zen = self._zenC()    # Type C = 'Doing X leads to Y'
            if mytype == 'D': zen = self._zenD()    # Type D = 'One X1's Y1 is an X2's Y2'
            ### ~ Tidy, join and return Zen as string ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while '' in zen: zen.remove('')
            ## Try to improve by removing repetition ##
            zen = string.join(zen)
            if mytype not in ['D']:
                wordlist = []
                for word in string.split(zen):
                    if len(word) < 5: continue
                    if word[:5].lower() in wordlist: return self.wisdom()
                    else: wordlist.append(word[:5].lower())
            for vowel in 'aeiou':
                zen = string.replace(zen,' a %s' % vowel,' an %s' % vowel)
                if zen.find('A %s' % vowel) == 0: zen = 'An %s' % zen[2:]
            return string.join(string.split(zen))
        except: self.log.errorLog('Bad vibes from Zen.wisdom()')
#########################################################################################################################
    ### <3> ### Random Wisdom Methods                                                                                   #
#########################################################################################################################
    def _zenA(self):    ### Generates Zen of the type 'The WISE MAN X BUT THE WWW XXX Y'
        '''Generates Zen of the type "The WISE MAN X BUT THE WWW XXX Y"'''
        if random.random() > 0.2: zen = [rje.randomList(['A','A','A','The','The','The','The','The'])[0], rje.randomList([self._adjective(),self._adverb()])[0], self._noun(), self._verb('A')]
        else: zen = ['The wise', self._noun(), self._verb('A')]
        if random.random() > 0.3: zen += [self._linker('A'), rje.randomList(['A','The'])[0].lower(), rje.randomList([self._adjective(),self._adverb()])[0], self._noun(), self._verb('A')]
        return zen
#########################################################################################################################
    def _zenB(self):    ### Generates Zen of the type 'It is better to X when Y'
        '''Generates Zen of the type "It is better X when Y"'''
        zen =  ['It is', self._adjective('B'), self._verb('B')]
        if random.random() > 0.5: zen += [self._linker('B'), rje.randomList(['A','The'])[0].lower(), rje.randomList([self._adjective(),self._adverb()])[0], self._noun(), self._verb('A')]
        return zen 
#########################################################################################################################
    def _zenC(self):    ### Generates Zen of the type 'Doing X leads to Y'
        '''Generates Zen of the type "Doing X leads to Y"'''
        zen = [self._verb('C'), rje.randomList(['nothing','anything',self._noun('Z')])[0],
               self._linker('C'), self._noun('C')]
        return zen
        #,rje.randomList(['A','The'])[0].lower()]
        #end = [[rje.randomList([self._adjective('C'),''])[0], self._noun('C')]] #,self._zenA()[-3:]]
        #return zen + rje.randomList(end)[0]
#########################################################################################################################
    def _zenD(self):    ### Generates Zen of the type 'One man's X is another man's Y'
        '''One man's X is another man's Y.'''
        zen = ['One',rje.randomList([self._adjective(),''])[0], self._noun()]
        zen[-1] = '%s\'s' % zen[-1]
        zen.append(self._noun(rje.randomList(['A','C','C','of','of'])[0]))
        if random.random() < 0.99: zen += ['is another',zen[2]]
        else: zen += ['is','%s\'s' % self._noun('Z')]
        zen.append(self._noun(rje.randomList(['A','C','C','of','of'])[0]))
        return zen
#########################################################################################################################
    ### <4> ### Zen Component Methods                                                                                   #
#########################################################################################################################
    def _adjective(self,ztype='A'): ### Returns a random adjective
        '''Returns a random adjective.'''
        alist = ['altruistic','benchmarking', 'better', 'bold', 'brave', 'cheesy', 'cocky', 'cowardly', 'crazy','coding',
                 'deadly', 'deceitful', 'degenerate',
                 'delicate', 'disillusioned', 'dynamic', 'existential','extraordinary',
                 'foolhardy', 'foolish','frantic','freaky','greedy','honest',
                 'kind','monastic','old','ordinary',
                 'passionate', 'pointless','prudent','quirky','quintessential',
                 'random', 'religious','RESTful','risky','running',
                 'scary','shrill','slippery', 'smooth', 'spiritual','semantic',
                 'stochastic','story-telling', 'superfluous', 'tasteless','unusual','usual',
                 'weird', 'wise', 'wise','young',
                 'zealous'] + ['zen'] * 3
        if ztype == 'A': alist += ['aerodynamic','aromatic','black','blue','Belgian','blessed',
                                   'cantankerous','cancerous','cursed','drug-addled','elongated',
                                   'Greek','Herculian','illiterate','killer','kindred',
                                   'opaque','orange','pale','pink','polymorphic','quality',
                                   'shell-like','short','strong','stringy','spongy','secular',
                                   'talented','tall','tiny','toxic','venomous',
                                   'white','yellow','youthful','zombie']
        if ztype == 'B': alist += ['folly','worse','nonsense','death','risky']
        return rje.randomList(alist)[0]
#########################################################################################################################
    def _noun(self,ztype='A'):    ### Returns a random noun
        '''Returns a random noun.'''
        nlist = ['ant', 'atheist', 'athelete','armadillo','assassin',
                 'badger', 'beetle', 'Bishop', 'boy', 'bushbaby', 'butterfly','cheerleader','billionaire','bilby',
                 'cat','cheese', 'chemist','child', 'cloud', 'communist', 'computer scientist','Creationist','Donald',
                 'diplomat', 'doctor', 'dragon', 'duck','elf', 'firefly',
                 'fish', 'fool', 'freak','fruit','fruitfly', 'fungus',
                 'girl', 'heretic','hound', 'jellyfish', 'kangaroo', 'knight', 'lady', 'ladybird',
                 'ladyboy','lion','kitten','kingfisher','king','killer',
                 'journalist',
                 'marsupial', 'mind', 'monkey', 'monster', 'mosquito', 'misogynist', 'nazi',
                 'philosopher', 'pig','pirate','postgrad',
                 'President','Professor','priest', 'prince', 'princess', 'programmer','python',
                 'queen', 'rabbit', 'Republican','runner',
                 'samurai', 'scientist', 'shrew', 'slug', 'snail',
                 'snake', 'soldier', 'spider', 'student', 'syntax error','smuggler',
                 'teapot', 'teenager', 'termite', 'theologian','theoretical physicist', 'Thespian','tiger','tiger snake',
                 'toddler',
                 'wallaby',
                 'warrior', 'wren','yeast','zealot','zombie','zulu'] + ['man'] * 8 + ['woman'] * 5
        if ztype == 'C':
            nlist = ['misery','happiness','poverty','wisdom','enlightenment','zen','dreams','passion','lunacy','plenty',
                     'alternative facts',
                     'death','disease','discovery','Chaos','religion','peril','philosophy','the Soul','debuggery','youth']
        if ztype == 'Z': return string.join([rje.randomList(['A','The'])[0].lower(),rje.randomList([self._adjective('A'),self._adverb('A'),''])[0], self._noun('A')])
        if ztype == 'of' or (ztype in ['A','B'] and random.random() < 0.2):
            nlist = nlist + ['shroom','pie','wine','gravy','egg','chocolate','cheese','banana','stool','horn','custard',
                             'tea','teacup','mug','motif', 'gene','genome', 'bucket','bucket','bucket',
                             'kebab','fruit','falafel','teapot',
                             'coffee','river','cake','cookie','sponge','abundance','repository','collection','library'] * 2
            nlist += ['cgi-bin']
            return string.join([rje.randomList(nlist)[0],'of',self._noun('C')])
        return rje.randomList(nlist)[0]
#########################################################################################################################
    def _linker(self,ztype='A'):    ### Returns a random linker
        '''Returns a random linker.'''
        zlist = []
        if ztype in ['A','B']: zlist += ['while','if','unless','when','because','whether or not']
        if ztype == 'A': zlist += ['but']
        #x#if ztype == 'B': zlist
        if ztype == 'C': zlist += ['leads to','leads to','shows','bamboozles','deceives',
                                   'enlightens','enlightens','enriches','enriches','exemplifies',
                                   'compresses','guzzumps','scoops',
                                   'invigorates','invites','rains on','shatters','tweaks','destroys','disturbs',
                                   'promotes','chastens','mocks','rejects','surprises','tests','unlocks','unravels',
                                   'juggles']
        return rje.randomList(zlist)[0]
#########################################################################################################################
    def _adverb(self,ztype='A'):    ### Returns a random adverb
        '''Returns a random linker.'''
        zlist = ['running','thinking','whistling','procrastinating','spanking','exhausted','spinning','talking',
                 'buzzing','drugged','skiing']
        return rje.randomList(zlist)[0]
#########################################################################################################################
    def _verb(self,ztype='A'):    ### Returns a random verb
        '''Returns a random verb.'''
        zlist = []
        znoun = self._noun('Z')
        if ztype == 'A': zlist += ['ponders','fishes','blows bubbles','scratches','fiddles','spends money',
                                   'cultivates %s' % self._noun('C'),'procrastinates','programs','poops',
                                   'consumes a %s' % self._noun('of'),'bathes in %s' % self._noun('C'),
                                   'cultivates %s' % self._noun('C'),'ponders %s' % self._noun('C'),
                                   'fiddles with %s' % self._noun('C')]
        if ztype == 'A': zlist += ['%s %s' % (rje.randomList(['marvels at','bitchslaps','tickles','worships','cultivates'])[0],znoun)]
        if ztype == 'B':
            zlist = ['to press flowers','to hop on one leg','to pee into the wind','to suck on a lollipop','to dream',
                     'to love','to procrastinate','to gesticulate','to Wii','to spout %s' % self._noun('C'),
                     'to shatter %s' % self._noun('C'),'to play jenga with blocks of %s' % self._noun('C'),
                     'to frolic through fields of %s' % self._noun('C'),
                     'to eat bananas','to program','to get your freak on','to high-five a %s' % self._noun(),
                     'to bathe in %s' % self._noun('C'),'to poke %s' % self._noun('Z'),
                     'to sequence a %s genome' % self._noun(),
                     'to sample the %s' % self._noun('of'),
                     'to drink from the %s' % self._noun('of'),
                     'to attend a conference about %s' % self._noun('C')]
            for i in range(3): zlist += ['to %s %s' % (rje.randomList(['smile at','shove','squeeze','punch','slap','stroke','tickle','shake','cuddle','worship','cultivate','nibble','defenestrate','curse'])[0],znoun)]
        if ztype == 'C': zlist += ['Doing','Eating','Slapping','Spanking','Loving','Poking','Stroking']
        return rje.randomList(zlist)[0]
#########################################################################################################################
### End of SECTION II: Zen Class                                                                                        #
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
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
        
    ### Rest of Functionality... ###
    try: Zen(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
