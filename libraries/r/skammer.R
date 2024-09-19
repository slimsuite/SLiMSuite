########################################################
### Skammer self-kmer assembly masker          ~~~~~ ###
### VERSION: 0.1.0                             ~~~~~ ###
### LAST EDIT: 20/08/24                        ~~~~~ ###
### AUTHORS: Richard Edwards 2024              ~~~~~ ###
### CONTACT: rich.edwards@uwa.edu.au           ~~~~~ ###
########################################################

# Self-Kmer Assembly Masker

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial version based on template.R
program <- "Skammer"
version <- "v0.1.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript skammer.R seqin=FILE [kfile=FILE] [basefile=STR] [threads=INT] [kmask=INT] [softmask=T/F] [unmask=T/F]
# : seqin=FILE = Genome assembly in fasta format
# : kfile=FILE = Self-kmer table (one value per base, pseudo-fasta format)
# : basefile=STR = Prefix of output files. [skammer]
# : threads=INT = Number of threads to run KAT with [8]
# : kmask=INT = Mask self-kmer of at least INT [4]
# : softmask=T/F = Whether to soft mask the assembly (lower case) rather than hard mask (N) [False]
# : unmask=T/F = Whether to unmask the assembly to upper case prior to masking [False]
# : smoothwin=INT = Smooth kmer values +/- INT bp [50]

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# os.popen('Rscript {0}skammer.R {1}'.format(rdir, optionstr)).readlines()

####################################### ::: OUTPUTS ::: ############################################
#!# List the outputs here

####################################### ::: TODO ::: ############################################
#!# Add parallel processing of regions if speedup needed.

####################################### ::: SETUP ::: ############################################
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(seqin="",kfile="",basefile="skammer",
                kmask=4,smoothwin=50,rdir="",threads=8,
                softmask=FALSE,unmask=FALSE,
                debug=FALSE,dev=FALSE,rscript=TRUE,
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
for(cmd in c("kmask","smoothwin","threads")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
for(cmd in c()){
  settings[[cmd]] = as.numeric(settings[[cmd]])
}
# for(cmd in c("listcmd")){
#   if(length(settings[[cmd]]) > 0){
#     settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
#   }else{
#     settings[[cmd]] <- settings[[cmd]]
#   }
# }
for(cmd in c("debug","dev","softmask","unmask")){
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

####################################### ::: FUNCTIONS ::: ############################################

##### ======================== Loading data functions ======================== #####

### ~ Load Depth File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# seqVector(filename,intvec=TRUE) - Returns list of seqname=vector
#i# loadTable(filename) 

################################## ::: RUN CODE ::: #######################################

##### ======================== Report key inputs ======================== #####

### ~ Key inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Assembly
seqin = settings$seqin
if(file.exists(seqin)){
  logWrite(paste("Assembly fasta file:",seqin))
}else{
  logWrite(paste("Assembly fasta file:",seqin,"not found!"))
  quit("no",2)  
}

#i# Kmer file of self-kmer counts per base
if(settings$kfile == ""){
  settings$kfile <- paste0(settings$basefile,".selfkat-counts.cvg")
}
kfile <- settings$kfile
logWrite(paste("Kmer File:",kfile))
# >> Make file with KAT if needed
if(! file.exists(kfile)){
  kbase <- paste0(settings$basefile,".selfkat")
  katcall <- paste("kat sect -t",settings$threads,"-o",kbase,seqin,seqin)
  system(katcall)
  if(kfile != paste0(kbase,"-counts.cvg")){
    success <- file.rename(paste0(kbase,"-counts.cvg"), kfile)
    if(success){
      logWrite(paste(paste0(kbase,"-counts.cvg"),"renamed to",kfile))
    }else{
      logWrite('#ERR Something went wrong!')
      quit("no",2)
    }
  }
}
# Check and quit if missing
if(! file.exists(kfile)){
  logWrite(paste("Kmer File:",kfile,"not found!"))
  quit("no",2)  
}
logWrite('#RCODE Setup complete.')


##### ======================== Load data ======================== #####

### ~ Load Kmer data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# 1. Load the self kmer list
klist = seqVector(kfile,shortname=TRUE)
settings$seqnames = names(klist)

### ~ Load Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# 2. Load the genome assembly
seqlist = fastaToList(seqin) 

#i# compare settings$seqnames to names(seqlist)
seqnames <- settings$seqnames[settings$seqnames %in% names(seqlist)]
logWrite(paste("#SEQIN",length(seqnames),"sequences loaded with self-kmer data."))

         
##### ======================== Process kmers ======================== #####

### ~ Smooth kmer data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(settings$smoothwin > 0){
  settings$smoothwin <- 0
  logWrite("Smoothing not yet implemented.")
}

if(settings$smoothwin > 0){
  
seqname <- seqnames[1]
kmers <- klist[[seqname]]

kwin <- settings$smoothwin
kL <- length(kmers)
kn <- (kwin * 2) + 1

kV <- c(rep(kmers[1],kwin),kmers,rep(kmers[kL],kwin))
cat(":", file = stderr())
for(i in 1:kwin){
  #logWrite(paste(i,kL+i-1,"+",kwin+i+1,kL+kwin+i))
  cat(".", file = stderr())
  kmers <- kmers + kV[i:(kL+i-1)] + kV[(kwin+i+1):(kL+kwin+i)]
}
kmers <- kmers / kn

cat("\n", file = stderr())
}

##### ======================== Output fasta ======================== #####

### ~ Add masking  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
maskbp <- 0
fullbp <- 0
#!# Add parallel processing
cat("Masking sequences", file = stderr())
seqdesc <- list()
for(seqname in seqnames){
  seqvec <- unlist(strsplit(seqlist[[seqname]], split=""))
  seqlen <- length(seqvec)
  if(settings$unmask){
    seqvec <- toupper(seqvec)
  }
  # Buffer the kmer vector
  kmers <- klist[[seqname]]
  klen <- length(kmers)
  kdiff <- seqlen - klen
  kstart <- as.integer(kdiff/2)
  kend <- kdiff - kstart
  kmers <- c(rep(kmers[1],kstart), kmers, rep(kmers[klen],kend))
  # Mask 
  if(settings$softmask){
    seqvec[kmers >= settings$kmask] <- tolower(seqvec[kmers >= settings$kmask])
  }else{
    seqvec[kmers >= settings$kmask] <- "N"
  }
  seqlist[[seqname]] <- paste0(seqvec,sep="",collapse="")
  kmaskn <- sum(kmers >= settings$kmask)
  maskbp <- maskbp + kmaskn
  fullbp <- fullbp + seqlen
  seqdesc[[seqname]] <- paste("Masked",kmaskn,"bp with kmer frequency >=",settings$kmask,paste0("(",signif(100*kmaskn/seqlen,3),"%)"))
  cat(".", file = stderr())
}
cat("\n", file = stderr())
masktxt <- "masked"
if(settings$softmask){
  masktxt <- "soft-masked"
}
if(settings$softmask){
  masktxt <- paste("pure",masktxt)
}  

logWrite(paste("#MASK Sequences",masktxt,"for kmer freq >=",settings$kmask,":",maskbp,"bp",paste0("(",signif(100*maskbp/fullbp,3),"%)"),"masked."))

### ~ Save fasta  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
filename <- paste0(settings$basefile,".skammer.fasta")
listToFasta(filename,seqlist,seqdesc,append=FALSE)
logWrite(paste("#SAVE",length(seqlist),masktxt,"sequences output to",filename))



### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE Skammer.R finished.")
if(settings$rscript){
  quit("no",0)
}
