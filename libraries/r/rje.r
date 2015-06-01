############# ::: GENERAL SETUP ::: #####################################
#library(Cairo)
#rdir = "/rhome/re1u06/Serpentry/"
#rdir = "/home/redwards/Serpentry/"
#rdir = "/data/ben/Serpentry/"
#rdir = "/scratch/RJE_Filestore/SBSBINF/Tools/_Serpentry/"
#rdir = "/home/re1u06/researchfiles/SBSBINF/Tools/_Serpentry/"
#rdir = "/home/re1u06/researchfiles/SBSBINF/Tools/svn/libraries/r/"
rdir = "/Users/redwards/Dropbox/_Repository_/slimsuite/libraries/r/"
argsettings = list("rdir"=rdir)

## ~ Read commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
args <- commandArgs(TRUE)
rtype = args[1]             #!# Update to use CommandArg instead?
arglist = list(rtype=rtype,basefile=args[2])    # Access with arglist[["basefile"]] or arglist$basefile
argsplit = strsplit(args[-c(1:2)],"=")
if(length(argsplit)>0){
  for(i in 1:length(argsplit)){
    cmd = argsplit[[i]][1]
    if(length(argsplit[[i]]) == 1){val = "T"}
    else{ val = argsplit[[i]][2] }
    arglist[[cmd]] = val
  }
}

## ~ General functions and objects~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
rdir = argsettings$rdir             # Can now be replaced by commandArgs
rjesource = function(rfile){
    source(paste(rdir,rfile,sep=""))
}
rjesource("rje_col.r")
rjesource("rje_misc.r")

############# ::: PARAMETER SETUP ::: ###################################


## ~ Commandline formatting functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
argBool = function(args,bools){   # Sets all values to TRUE/FALSE
  for(cmd in bools){
    if(sum(c('f','F') == (substr(args[[cmd]],1,1))) > 0){ args[[cmd]] = FALSE }
    else{ args[[cmd]] = TRUE }
  }
  return(args)
}
argNum = function(args,nums){   # Converts all values to numerics
  for(cmd in nums){ args[[cmd]] = as.numeric(args[[cmd]]) }
  return(args)
}

############# ::: CALL APPROPRIATE SCRIPT ::: ###########################
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
