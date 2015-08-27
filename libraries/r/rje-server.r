############# ::: GENERAL SETUP ::: #####################################
library(Cairo)
rdir = "/home/re1u06/Tools/svn/libraries/r/"
#rdir = "/home/redwards/Serpentry/"
#rdir = "/data/ben/Serpentry/"
#rdir = "/scratch/RJE_Filestore/SBSBINF/Tools/_Serpentry/"
#rdir = "/home/re1u06/researchfiles/SBSBINF/Tools/_Serpentry/"
#rdir = "/home/re1u06/researchfiles/SBSBINF/Tools/svn/libraries/r/"
#rdir = "/Users/redwards/Dropbox/_Repository_/slimsuite/libraries/r/"
rjesource = function(rfile){
    source(paste(rdir,rfile,sep=""))
}
rjesource("rje_col.r")
rjesource("rje_misc.r")

############# ::: PARAMETER SETUP ::: ###################################

## ~ Read commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
args <- commandArgs(TRUE)
rtype = args[1]

############# ::: CALL APPROPRIATE SCRIPT ::: ###########################
slimjims = c("spokealn","motif")
for (sj in slimjims){
	if (rtype == sj){ rjesource("slimjim.r") }
}

if (rtype == "interactome" | rtype == "dispng" | rtype == "nested" | rtype == "pureppi"){ rjesource("rje_ppi.r") }

if (rtype == "hminteractome" | rtype == "hmnested"){ rjesource("hmhtml.r") }

if (substr(rtype,1,2) == "go"){ rjesource("aphid_go.r") }

if (rtype == "sfmap2png") { rjesource("sfmap2png.r") }
if (rtype == "rel2png") { rjesource("sfmap2png.r") }
if (rtype == "occaln") { rjesource("sfmap2png.r") }

if (rtype == "sfslim2png") { rjesource("sfslim2png.r") }
