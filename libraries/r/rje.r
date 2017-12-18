########################################################
### RJE R Script controller                          ###
### VERSION: 0.1.0                             ~~~~~ ###
### LAST EDIT: 30/06/17                        ~~~~~ ###
### AUTHORS: Richard Edwards 2016              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

################# ::: HISTORY ::: ######################
# v0.1.0 : Added code for self-detecting rdir. Added History and header.

############# ::: GENERAL SETUP ::: #####################################
#library(Cairo)
#rdir = "/rhome/re1u06/Serpentry/"
#rdir = "/home/redwards/Serpentry/"
#rdir = "/data/ben/Serpentry/"
#rdir = "/scratch/RJE_Filestore/SBSBINF/Tools/_Serpentry/"
#rdir = "/home/re1u06/researchfiles/SBSBINF/Tools/_Serpentry/"
#rdir = "/home/re1u06/researchfiles/SBSBINF/Tools/svn/libraries/r/"
#rdir = "/Users/redwards/Dropbox/_Repository_/slimsuite/libraries/r/"
rdir = getSrcDirectory(function(dummy) {dummy})
#i# This gives the directory of the file where the statement was placed (where the dummy function is defined). It can then be used to set the working direcory and use relative paths e.g.
#i# Previously, rdir should be be read from the commandline as an option. This *should* remove that need.

#!# Add some if(interactive()){} toggles where it should have different behaviour if run non-interactively.

## ~ Read commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
args <- commandArgs(TRUE)
if(length(commandArgs(TRUE)) > 0){
  arglist = list(rdir=rdir,rtype=args[1],basefile=args[2])    # Access with arglist[["basefile"]] or arglist$basefile
}
rtype = arglist$rtype             
basefile = arglist$basefile
argsplit = strsplit(args[-c(1:2)],"=")
if(length(argsplit)>0){
  for(i in 1:length(argsplit)){
    cmd = argsplit[[i]][1]
    if(length(argsplit[[i]]) == 1){val = "T"}
    else{ val = argsplit[[i]][2] }
    arglist[[cmd]] = val
  }
}
writeLines("Commandline arguments:")
print(arglist)

## ~ General functions and objects~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
rdir = arglist$rdir    # Hopefully will have been set by commandline and over-ruled default.
rjesource = function(rfile){
    source(paste(rdir,rfile,sep=""))
}
rjesource("rje_col.r")
rjesource("rje_misc.r")

############# ::: PARAMETER SETUP ::: ###################################
## ~ Commandline formatting functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
argBool = function(args,bools){   # Sets all values to TRUE/FALSE
  for(cmd in bools){
    if(length(args[[cmd]])>0){
      if(sum(c('f','F') == (substr(args[[cmd]],1,1))) > 0){ args[[cmd]] = FALSE }
      else{ args[[cmd]] = TRUE }
    }
  }
  return(args)
}
argNum = function(args,nums){   # Converts all values to numerics
  for(cmd in nums){ 
    if(length(args[[cmd]])>0){
      args[[cmd]] = as.numeric(args[[cmd]]) 
    }
  }
  return(args)
}
argList = function(args,lists){
  for(cmd in lists){ 
    if(length(args[[cmd]])>0){
      args[[cmd]] = strsplit(args[[cmd]],",")[[1]] 
    }
  }
  return(args)
}
writeLines("rje.R setup complete")
writeLines("")

############# ::: CALL APPROPRIATE SCRIPT ::: ###########################
# Newer scripts that make use of commandline argument capability
if (rtype == 'pagsat'){ rjesource("pagsat.R") }
if (rtype == 'pagsat_V1'){ rjesource("pagsat_V1.R") }
if (rtype == 'snapper'){ rjesource("snapper.R") }
if (rtype == 'snpfreq'){ rjesource("snpfreq.R") }
if (rtype == 'pacbio'){ rjesource("pacbio.R") }
if (rtype == 'tree'){ rjesource("tree.R") }
if (rtype == 'samtools'){ rjesource("samtools.R") }
if (rtype == 'samdepth'){ rjesource("samtools.R") }
if (rtype == 'samreadlen'){ rjesource("samtools.R") }
if (rtype == 'samunzip'){ rjesource("samtools.R") }




# Old scripts, not yet implementing commands etc. properly
slimjims = c("spokealn","motif")
for (sj in slimjims){
	if (rtype == sj){ rjesource("slimjim.r") }
}

if (rtype == 'revert'){ rjesource("revert.r") }

if (rtype == "interactome" | rtype == "dispng" | rtype == "nested" | rtype == "pureppi"){ rjesource("rje_ppi.r") }

if (rtype == "hminteractome" | rtype == "hmnested"){ rjesource("hmhtml.r") }

if (substr(rtype,1,2) == "go"){ rjesource("aphid_go.r") }

if (rtype == "sfmap2png") { rjesource("sfmap2png.r") }
if (rtype == "rel2png") { rjesource("sfmap2png.r") }
if (rtype == "occaln") { rjesource("sfmap2png.r") }

if (rtype == "sfslim2png") { rjesource("sfslim2png.r") }

writeLines(c("|--------|","Warnings:"))
warnings()

