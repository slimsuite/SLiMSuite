########################################################
### CentroFish Centromere Repeat Fishing       ~~~~~ ###
### VERSION: 0.1.0                             ~~~~~ ###
### LAST EDIT: 20/08/24                        ~~~~~ ###
### AUTHORS: Richard Edwards 2024              ~~~~~ ###
### CONTACT: rich.edwards@uwa.edu.au           ~~~~~ ###
########################################################

# This script is for homology-based chromosome synteny plots

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial version based on OG14 data exploration
version = "v0.1.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript centrofish.R trashin=FILE [peakin=FILE] [replen=LIST] [basefile=FILE] [minwidth=INT] [maxgap=INT] [chrprefix=X] [optimise=T/F]

# : trashin=FILE = Summary file from TRASH [Summary.of.repetitive.regions.csv]
# : peakin=FILE = Optional file of TRASH peaks to pick out kmer lengths [peaks.csv]
# : replen=LIST = Optional list of repeat lengths to use as candidate CEN repeats []
# : basefile=FILE = Prefix for outputs [centrofish]
# : minwidth=INT = minimum width of a repeat block to use (bp) [0]
# : maxgap=INT = maximum gap between repeat blocks for merging (bp) [10000]
# : chrprefix=X = text prefix for chromosome-level scaffolds [SUPER]
# : optimise=T/F = whether to perform additional optimisation of minwidth [TRUE]
# : debug=T/F = whether to switch on additional debugging outputs [FALSE]
# : dev=T/F = whether to switch on dev mode during code updates [FALSE]

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# os.popen('Rscript {0}chromsyn.R {1}'.format(rdir, optionstr)).readlines()

#i# Usage within R:
# Set an override vector of commandline arguments: this will replace argvec read from the commandline
# Use source() to run the script:
# source("$PATH/chromsyn.R")

####################################### ::: OUTPUTS ::: ############################################
#!# List the outputs here

####################################### ::: TODO ::: ############################################
#!# Initial working version.
#!# 

####################################### ::: SETUP ::: ############################################
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(trashin="Summary.of.repetitive.regions.csv",
                peakin="peaks.csv", replen="",
                basefile="centrofish",
                minwidth=0,maxgap=10000,
                chrprefix="SUPER",optimise=TRUE,
                debug=FALSE,dev=FALSE,
                rdir="",runpath="",
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
    settings[[cmdv[1]]] = ""    
  }
}
#i# integer parameters
for(cmd in c("minwidth","maxgap")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
#i# other numeric parameter
for(cmd in c()){
  settings[[cmd]] = as.numeric(settings[[cmd]])
}
#i# list parameters
for(cmd in c("replen")){
  if(sum(grep(",",settings[[cmd]],fixed=TRUE)) > 0){
    settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
  }else{
    settings[[cmd]] <- settings[[cmd]]
  }
}
#i# logical parameters
for(cmd in c("debug","dev","optimise")){
  settings[[cmd]] = as.logical(settings[[cmd]])
}

#i# Set warnings based on debug
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
if(! settings$runpath == ""){
  setwd(settings$runpath)
}
logWrite(paste("#RCODE ChromSyn.R:",version))
logWrite(paste("#PATH Running from:",getwd()))
for(cmd in names(settings)[order(names(settings))]){
  logWrite(paste("CMD:",cmd,"=",paste0(settings[[cmd]],collapse=",")))
}
#dir.create(settings$plotdir, showWarnings = FALSE)

### ~ Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(! "tidyverse" %in% installed.packages()[,"Package"]){
  install.packages("tidyverse")
}
library(tidyverse)
library(RColorBrewer)
library(gtools)
settings$ggstatsplot = "ggstatsplot" %in% installed.packages()[,"Package"]
if(settings$ggstatsplot){
  library(ggstatsplot)
}
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


################################## ::: PLOTTING FUNCTIONS ::: #######################################


################################## ::: RUN CODE ::: #######################################

##### ======================== Report key inputs ======================== #####

### ~ Key inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# TRASH Summary File = vital input
sumfile <- settings$trashin
logWrite(paste("TRASH Summary File:",sumfile))
if(! file.exists(sumfile)){
  msg <- paste("Cannot find TRASH Summary file:",sumfile)
  if(settings$rscript){
    logWrite(msg)
    quit("no",2)  
  }else{
    stop(msg)
  }
}
#i# Peaks file
peakfile <- settings$peakin
if(! file.exists(peakfile)){
  msg <- paste("Cannot find TRASH Peaks file:",peakfile)
  peakfile <- NA
}
logWrite(paste("TRASH Peaks File:",peakfile))

logWrite('#RCODE Setup complete.')

##### ======================== Load data ======================== #####

### ~ Summary File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
sumdb <- read_csv(sumfile,show_col_types = FALSE) %>%
  rename(replen = most.freq.value.N, repnum = repeats.identified, consnum = consensus.count, repseq = consensus.primary) %>%
  select(name, start, end, width, replen, consnum, repnum, repseq) %>%
  arrange(name,start,end)
logWrite(paste('#LOAD Loaded',nrow(sumdb),"repeat blocks from",sumfile))

### ~ Peaks file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
peakdb <- tibble(repeat_length = c())
if(! is.na(peakfile)){
  peakdb <- read_csv(peakfile,show_col_types = FALSE)
}
logWrite(paste('#PEAKS',nrow(peakdb),"peaks loaded."))

if(settings$debug | settings$dev){
  head(peakdb)
  head(numdb)
}

##### ======================== Compile data ======================== #####

#i# Minimum width for a repeat to be considered for centromeres
minwidth <- settings$minwidth
#i# Minimum gap between repeats to join into longer repeat unit
maxgap <- settings$maxgap
#i# Prefix for chromosomes
chrprefix <- settings$chrprefix
#i# Output prefix
basefile <- settings$basefile

# Sort by replen, name, start
repdb <- sumdb %>% arrange(replen, name, start) %>%
  mutate(repid = paste(name,replen), chrom = str_starts(name, chrprefix))

# Filter to candidates
if(nrow(peakdb) > 0){
  logWrite(paste('Filter to peaks:',paste(peakdb$repeat_length,collapse=", ")))
  repdb <- repdb[repdb$replen %in% peakdb$repeat_length,]
  logWrite(paste('#PEAKS',nrow(repdb),"repeat blocks after filtering to",nrow(peakdb),"peaks."))
}
if(length(settings$replen) == 1 & settings$replen[1] == ""){
  settings$replen <- c()
}
logWrite(settings$replen)
if(length(settings$replen) > 0){
  logWrite(paste('Filter to repeat unit lengths:',paste(settings$replen,collapse=", ")))
  repdb <- repdb[repdb$replen %in% as.integer(settings$replen),]
  logWrite(paste('#REPLEN',nrow(repdb),"repeat blocks after filtering to",length(settings$replen),"repeat sizes."))
}
#i# Summarise
cat("Chromosome repeats:\n")
table(repdb$chrom)

##### ======================== Process repeats ======================== #####

### ~ Merge repeat blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# Combine lengths within distance
logWrite(paste('#MERGE Merging repeat blocks upto',maxgap,"bp apart."))
repdb$merged <- FALSE
merging <- TRUE
while(merging){
  merging <- FALSE
  for(j in 2:nrow(repdb)){
    i <- j - 1
    # Check for merging rules
    if(repdb$merged[i]){ next }
    if(repdb$repid[i] != repdb$repid[j]){ next }
    if((repdb$end[i] + maxgap) < repdb$start[j]){ next }
    # Merge i into j
    merging <- TRUE
    repdb$merged[i] <- TRUE
    repdb$start[j] <- repdb$start[i]
    repdb$width[j] <- repdb$width[i] + repdb$width[j]
    repdb$consnum[j] <- repdb$consnum[i] + repdb$consnum[j]
    repdb$repnum[j] <- repdb$repnum[i] + repdb$repnum[j]
    
  }
}
logWrite(paste("#MERGE",sum(repdb$merged),"repeat blocks merged and removed."))
repdb <- repdb %>% filter(! merged)

### ~ Optimise min repeat block size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#?# Should we add a buffer for kmers of similar size?
chrnum <- nrow(unique(repdb %>% filter(chrom) %>% select(name)))
logWrite(paste("#CHROM",chrnum,"Chromosome scaffolds with retained repeat blocks"))
repdb <- repdb %>% arrange(width) %>% select(-repid, -merged) %>%
  filter(width >= minwidth)
if(settings$optimise & chrnum > 0){
  logWrite(paste('Optimising min width to retain blocks in', chrnum, 'chromosomes...'))
  # Establish a minimum length that maximises the number of chromosomes with a centromere
  for(rlen in unique(repdb$width)){
    if(nrow(unique(repdb %>% filter(chrom, width >= rlen) %>% select(name))) < chrnum){
      break
    }
    minwidth <- rlen
  }
}

### ~ Filter on min repeat block size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite(paste("#MINREP Min width for candidate centromere:",signif(minwidth/1000,4),"kb"))
repdb <- repdb %>% filter(width >= minwidth) %>% arrange(name,start)
logWrite(paste('#REPLEN',nrow(repdb),"repeat blocks after filtering to",signif(minwidth/1000,4),"kb width."))

if(settings$debug){
  table(repdb$replen)
  chromrep <- unique(repdb %>% filter(chrom) %>% select(name,replen))
  table(chromrep$replen)
}

# Identify scaffolds that are mostly centromere
#!# Need to load lengths for this


##### ======================== Generate Output ======================== #####

### ~ Summary Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
outfile = paste(basefile,"centrofish.tsv",sep=".",collapse=".")
logWrite(paste("#SAVE",nrow(repdb),"Filtered repeats output to",outfile))
write.table(repdb,outfile,sep="\t",quote=FALSE,row.names=FALSE)


### ~ ChromSyn Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# Output Shape=-2 for pairs of triangles +25 -24 to add to chromsyn plots
# SeqName,Start,End,Strand,Col,Shape
maxrep <- max(repdb$replen)
fwddb <- repdb %>% mutate(Strand="+",Col=if_else(replen==maxrep,"darkblue","darkgreen"),Shape=-2) %>%
  rename(SeqName = name, Start = start, End = end) %>%
  select(SeqName,Start,End,Strand,Col,Shape)
bwddb <- repdb %>% mutate(Strand="-",Col=if_else(replen==maxrep,"darkblue","darkgreen"),Shape=-2) %>%
  rename(SeqName = name, Start = start, End = end) %>%
  select(SeqName,Start,End,Strand,Col,Shape)
cenrepdb <- bind_rows(fwddb,bwddb)
outfile = paste(basefile,"cenrep.csv",sep=".",collapse=".")
logWrite(paste("#SAVE",nrow(repdb),"Repeats for ChromSyn plotting output to",outfile))
write.table(cenrepdb,outfile,sep=",",quote=FALSE,row.names=FALSE)


##### ======================== Tidy and Finish ======================== #####
if(file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE CentroFish.R finished.")
#quit("no",0)

