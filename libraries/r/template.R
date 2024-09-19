########################################################
### Template                                   ~~~~~ ###
### VERSION: 0.0.0                             ~~~~~ ###
### LAST EDIT: 27/05/22                        ~~~~~ ###
### AUTHORS: Richard Edwards 2022              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# Replace with one-line description of code purpose here.

####################################### ::: HISTORY ::: ############################################
# v0.0.0 : Initial version based on template.R
program <- "Template"
version <- "v0.0.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript template.R seqin=FILE depfile=FILE basefile=STR
# : seqin=FILE = Genome assembly in fasta format
# : depfile=FILE = Full depth data table (one value per base, pseudo-fasta format)
# : basefile=STR = Full BUSCO table. If a number, will interpret as scdepth value.
# : minhit=NUM = Min percentage coverage to be kept as a Hit [40]
# : scdepth=NUM = Single copy read depth to use in place of BUSCO-derived depth
# : buscocn=T/F = Whether to use the BUSCO-derived CN setting if possible to calculate.
# : threads=NUM = Number of threads to use. [0 will use detectCores() - 1] (Not yet implemented)

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# os.popen('Rscript {0}deppurgehap.R {1}'.format(rdir, optionstr)).readlines()

####################################### ::: OUTPUTS ::: ############################################
#!# List the outputs here

####################################### ::: TODO ::: ############################################
#!# Add parallel processing of regions if speedup needed.

####################################### ::: SETUP ::: ############################################
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(adjust=12,scdepth=0,busco="",threads=0,
                debug=FALSE,
                gablam='None',minhit=40,
                basefile="depthcopy",rscript=TRUE,
                outlog=stdout())
settings <- defaults
argvec = commandArgs(TRUE)
if("override" %in% ls()){
  argvec = override
  settings$rscript = FALSE
}
for(cmd in argvec){
  cmdv = strsplit(cmd,'=',TRUE)[[1]]
  if(length(cmdv) > 1){
    settings[[cmdv[1]]] = cmdv[2]
  }else{
    settings[[cmdv[1]]] = TRUE    
  }
}
for(cmd in c("adjust","threads")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
for(cmd in c("scdepth","minhit")){
  settings[[cmd]] = as.numeric(settings[[cmd]])
}
for(cmd in c("listcmd")){
  if(length(settings[[cmd]]) > 0){
    settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
  }else{
    settings[[cmd]] <- defaults[[cmd]]
  }
}
for(cmd in c("debug")){
  settings[[cmd]] = as.logical(settings[[cmd]])
}

oldwarn <- getOption("warn")

if(settings$debug){
  writeLines(argvec)
}else{
  options(warn = -1)
}

### ~ logWrite function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite <- function(logstr){
  writeLines(paste0("[",date(),"] ",logstr),con=settings$outlog)
}
logWrite(paste("#RCODE",program,version))
for(cmd in names(settings)[order(names(settings))]){
  if(cmd %in% c("fields")){
    for(x in names(settings[[cmd]])){
      logWrite(paste("CMD:",cmd,":",x,"=",paste0(settings[[cmd]][[x]],collapse=",")))
    }
  }else{
    logWrite(paste("CMD:",cmd,"=",paste0(settings[[cmd]],collapse=",")))
  }
}

### ~ Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(! "tidyverse" %in% installed.packages()[,"Package"]){
  install.packages("tidyverse")
}
library(tidyverse)
settings$writexl = "writexl" %in% installed.packages()[,"Package"]
if(settings$writexl){
  library(writexl)
}else{
  logWrite("#XLXS Install writexl package for compiled Excel file output.")
}

### ~ Load R libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
if(settings$rdir == ""){
  settings$rdir <- getScriptPath()
}
sfile <- paste0(settings$rdir,"/rje_load.R")
logWrite(sfile)
source(sfile)

### ~ Setup parallelisation and threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Currently no parallelisation, but remove this line when available
#i# Use clusterExport(settings$cl, "<env_object>") to get the data into the cluster
#i# USe clusterApply(settings$cl, <vector>, <function>) to run function in parallel
settings$threads <- 1
if(settings$threads != 1 & "parallel" %in% installed.packages()[,"Package"]){
  library(parallel)
  if(settings$threads == 0){
    settings$threads <- detectCores() - 1
  }
  if(settings$threads < 0){
    settings$threads <- 1
  }
  if(settings$threads > 1){
    settings$cl <- makeCluster(settings$threads)
    logWrite(paste("#PARA Local cluster with",settings$threads,"threads created."))
  }
}else{
  if(settings$threads != 1){
    logWrite("#PARA Install parallel package for multithreaded run.")
  }
  settings$threads <- 1
}

####################################### ::: FUNCTIONS ::: ############################################

##### ======================== Loading data functions ======================== #####


################################## ::: PLOTTING FUNCTIONS ::: #######################################


################################## ::: RUN CODE ::: #######################################

##### ======================== Report key inputs ======================== #####


##### ======================== Tidy and Finish ======================== #####
if(file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE ChromSyn.R finished.")