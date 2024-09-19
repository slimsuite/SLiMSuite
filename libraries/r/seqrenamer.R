########################################################
### SeqRenamer: Sequence identifier mapper     ~~~~~ ###
### VERSION: 0.3.0                             ~~~~~ ###
### LAST EDIT: 10/09/24                        ~~~~~ ###
### AUTHORS: Richard Edwards 2022              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for mapping/updating sequence identifiers

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial version based on SynBad data exploration
# v0.1.1 : Updated for different transformation table formats.
# v0.1.2 : Added option to remove edit from output. Fixed pos bug.
# v0.2.0 : 
# v0.3.0 : Added split=T/F option to split any regions divided by transformations.
version = "v0.3.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript seqrenamer.R transform=FILE in=FILE [seqlen=FILE] [out=FILE] [basefile=FILE] [fields=DICT]

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# os.popen('Rscript {0}seqrenamer.R {1}'.format(rdir, optionstr)).readlines()

#i# Usage within R:
# Set an override vector of commandline arguments: this will replace argvec read from the commandline
# Use source() to run the script:
# source("$PATH/seqrenamer.R")

####################################### ::: INPUTS ::: ############################################
# The primary inputs are a translation table [transform=FILE] and a table to be converted [in=FILE]

## transform=FILE
# > D$transdb : data frame with SeqName, [Start, End,] NewSeqName, [Shift, RevComp]
#i# If Start and End are missing, will be 1,-1 

# in=FILE
# > D$in : data frame with SeqName, [Pos, Start, End, Strand]

# fields=DICT : list of field:fieldname, separated by commas (fieldname is in input)
# > D$fields : list object

# seqlen=FILE : File of SeqName, SeqLen
# > D$seqlen : data frame of SeqName and SeqLen
#i# Needed for reverse complementing positions. Can be in the input table.

####################################### ::: OUTPUTS ::: ############################################
#!# List the outputs here

####################################### ::: TODO ::: ############################################
#!# Basic version that takes input and translation table and generates output.
#!# Add position shift during mapping.
#!# Add strand and sequence length to swap orientation.
#!# Add start and end swap for strand change.
#!# Options for inclusive and exclusive joins of inputs.
#?# Add field case toggle for pure lower case fields?

####################################### ::: SETUP ::: ############################################
# Rscript seqrenamer.R transform=FILE in=FILE [seqlen=FILE] [filter=T/F] [out=FILE] [basefile=FILE] [fields=DICT] [camelcase=T/F]
# : transform=FILE = Transformation table with SeqName, [Start, End,] NewSeqName, [Shift, RevComp]
# : in=FILE = Input file to be modified SeqName, [Pos, Start, End, Strand, SeqLen]
# : seqlen=FILE = File of SeqName and SeqLen, if SeqLen not in input and RevComp in transform=FILE
# : filter=T/F = Whether to only output sequences/regions in transform table [True]
# : split=T/F = Whether to split any regions divided by transformations [True]
# : out=FILE = Output file with modified sequence names
# : fields=DICT = Dictionary of field:fieldname, where fieldname is the column name in the input FILE
# : camelcase=T/F = Whether final output should be in CamelCase (True), or lower case (False) [True]
# : addedit=T/F = Whether to include the edit field in the output [True]
# : basefile=FILE = Prefix for outputs [seqrenamer]
# : debug=T/F = whether to switch on additional debugging outputs [FALSE]
# : dev=T/F = whether to switch on dev mode during code updates [FALSE]
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(transform="transform.tsv",`in`="in.tsv",out="out.tsv",seqlen="",
                rscript=TRUE,debug=FALSE,dev=FALSE,filter=TRUE,addedit=TRUE,
                basefile="seqrenamer",fields=list(),camelcase=TRUE,split=TRUE,
                rdir="",
                outlog=stdout())

settings <- defaults
argvec = commandArgs(TRUE)
#override <- c("transform=WatersSup2.txt","in=stats/Chicken.CHICK.500kb.sequences.tdt","out=CHICKCHR.sequences.tsv","fields=seqname:name")
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
for(cmd in c()){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
for(cmd in c()){
  settings[[cmd]] = as.numeric(settings[[cmd]])
}
for(cmd in c()){
  if(length(settings[[cmd]]) > 0){
    settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
  }else{
    settings[[cmd]] <- defaults[[cmd]]
  }
}
for(cmd in c("fields")){
  if(length(settings[[cmd]]) > 0){
    #i# Lowercase mapping
    cdict <- list()
    for(el in strsplit(settings[[cmd]],',',TRUE)[[1]]){
      argval <- tolower(strsplit(el,':',TRUE)[[1]])
      cdict[[argval[1]]] <- argval[[2]]
    }
    settings[[cmd]] <- cdict
  }else{
    settings[[cmd]] <- defaults[[cmd]]
  }
}
for(cmd in c("debug","dev","camelcase","filter","addedit","split")){
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
logWrite(paste("#RCODE SeqMapper.R:",paste(version,collapse=" ")))
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

### ~ Transform genomic locations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# > D : data frame with seqname, [start, end, pos] and strand
# > transdb : data frame with seqname, start, end, seqname, shift, revcomp
# > strict=FALSE : boolean as to whether to exclude any regions in D that are not in transdb
# > logwrite=TRUE : Whether to summarise to logWrite()
# [ ] : Future > split=TRUE : whether to split partially overlapping regions rather than dumping them
transformLocations <- function(D,transdb,strict=FALSE,logwrite=TRUE){
  # Set up fields to remove later  
  dumpfield <- c()
  if(! settings$addedit){
    dumpfield <- c("edit")
  }
  if(! "seqname" %in% colnames(D)){
    dumpfield <- c(dumpfield,"seqname")
    for(pfield in c("id","contig")){
      if(pfield %in% colnames(D)){
        D$seqname <- D[[pfield]]
      }
    }
  }
  if(! "seqname" %in% colnames(D)){
    logWrite("Failed to find seqname field!")
  }
  # Debugging
  if(settings$debug){
    logWrite(paste(colnames(D),collapse=" "))
    logWrite(paste(unique(D$seqname),collapse=" "))
  }
  # Check that unique(seqname, start, end) is unique
  transdb <- transdb %>% filter(! seqname == "Gap")
  ucheck <- transdb %>% select(seqname, start, end)
  ucheck <- unique(ucheck)
  if(nrow(ucheck) != nrow(transdb)){
    abort("transformLocations() was given transdb with non-unique seqname, start, end data")
  }
  # Add start, end, strand and/or pos for processing
  pfield <- "strand"
  if(! pfield %in% colnames(D)){
    dumpfield <- c(dumpfield,pfield)
    D[[pfield]] <- "+"
  }
  pfield <- "pos"
  if(! pfield %in% colnames(D)){
    dumpfield <- c(dumpfield,pfield)
    if("window" %in% colnames(D)){
      D[[pfield]] <- D$window
    }else{
      D[[pfield]] <- D$start
    }
  }
  for(pfield in c("start","end")){
    if(! pfield %in% colnames(D)){
      dumpfield <- c(dumpfield,pfield)
      D[[pfield]] <- D$pos
    }
  }
  # Make sure everything matches if no maxlen given
  maxlen <- max(D$end,D$start,D$pos,na.rm=TRUE)
  #if(max(transdb$end) < 1){
  #  transdb$end <- maxlen
  #}
  # Work through D and shift
  #i# D is the table of regions to be shifted
  #i# transdb is the table of shifts to apply
  newD <- D
  newD$edit <- FALSE
  for(i in 1:nrow(transdb)){
    reg <- transdb[i,]
    if(reg$seqname == reg$newseqname & reg$shift == 0 & ! reg$revcomp){ next }
    if(reg$end < 1){
      if("seqlen" %in% names(DF) & reg$seqname %in% DF$seqlen$SeqName){
        reg$end <- DF$seqlen[DF$seqlen$SeqName == reg$seqname,]$SeqLen
      }else{
        reg$end <- maxlen
      }
    }
    #i# Construct a vector of TRUE or FALSE for D rows that need shifting
    newDreg <- D$seqname == reg$seqname & ( (D$start >= reg$start & D$end <= reg$end) | (reg$start == 1 & is.na(D$start)))
    if(settings$debug){
      writeLines(paste(sum(newDreg),reg$seqname,"transformation(s) ->",reg$newseqname))
    }
    if(sum(newDreg) == 0){ next }
    if(sum(newD[newDreg,]$edit) > 0){
      abort("transformLocations() is trying to edit same region multiple times")
    }
    newD[newDreg,]$seqname <- reg$newseqname
    if(reg$revcomp){
      # Need to move everything to be from reg$end not 1 and swap
      newD[newDreg,]$start <- reg$end - D[newDreg,]$end + reg$start + reg$shift
      newD[newDreg,]$end <- reg$end - D[newDreg,]$start + reg$start + reg$shift
      newD[newDreg,]$pos <- reg$end - D[newDreg,]$pos + reg$start + reg$shift
      if("window" %in% colnames(D)){
        newD[newDreg,]$window <- reg$end - D[newDreg,]$window + reg$start + reg$shift
      }
      if(sum(newDreg & D$strand == "+")){ #   ("+" %in% newD[newDreg,]$strand){
        newD[newDreg & D$strand == "+",]$strand <- "-"
      }
      if(sum(newDreg & D$strand == "-")){ #if("-" %in% newD[newDreg,]$strand){
        newD[newDreg & D$strand == "-",]$strand <- "+"
      }
    }else{
      # Otherwise, simple shift
      newD[newDreg,]$start <- D[newDreg,]$start + reg$shift
      newD[newDreg,]$end <- D[newDreg,]$end + reg$shift
      newD[newDreg,]$pos <- D[newDreg,]$pos + reg$shift
      if("window" %in% colnames(D)){
        newD[newDreg,]$window <- D[newDreg,]$window + reg$shift
      }
    }
    newD[newDreg,]$edit <- TRUE
  }
  # Summarise
  #if(logwrite & "logWrite" %in% ls()){
  logWrite(paste(sum(newD$edit),"of",nrow(D),"regions mapped to new names/locations"))
  #}
  # Drop dumpfields
  if(settings$debug){
    logWrite(paste(colnames(D),collapse=" "))
    logWrite(paste(dumpfield,collapse=" "))
  }
  for(pfield in dumpfield){
    newD <- select(newD,-!!pfield)
  }
  logWrite('Transformation complete')
  return(newD)
}

####################################### ::: BACKBONE ::: ############################################
#1# Read input data and convert fields to lower case. Generate camel list of lower:CamelCase conversion
D <- list()  # Input data will be loaded here
#i# Conversion of lower case to camel case fields
camel <- list(seqname="SeqName",start="Start",end="End",pos="Pos",seqlen="SeqLen",strand="Strand")
defval <- list(start=1,end=-1,shift=0,revcomp=FALSE)


### ~ Load data ~~~~~~~~~~~~~~~~~~~~~~ ###

#i# Transformation table
## transform=FILE
# > D$transdb : data frame with SeqName, [Start, End,] NewSeqName, [Shift, RevComp]
#i# If Start and End are missing, will be 1,-1 
delimit <- "\t"
filename <- settings$transform
if(endsWith(filename,".txt")){
  delimit <- " "
}
if(endsWith(filename,".csv")){
  delimit <- ","
}
if(file.exists(filename)){
  D$transdb <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="")
  for(fname in colnames(D$transdb)){
    camel[[tolower(fname)]] <- fname
  }
  colnames(D$transdb) <- tolower(colnames(D$transdb))
  for(field in names(defval)){
    if(! field %in% colnames(D$transdb)){
      D$transdb[[field]] <- defval[[field]]
    }
  }
  logWrite(paste('#LOAD',nrow(D$transdb),"transformations loaded from",filename))
}else{
  if(settings$rscript){
    quit("no",2)  
  }else{
    stop(paste("Cannot find transformation table:",filename))
  }
}

#i# Input table
## in=FILE
# > D$input : data frame with SeqName, [Pos, Start, End, Strand]
filename <- settings$`in`
if(file.exists(filename)){
  #D$input <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="")
  TD <- loadTable(filename)
  D$input <- TD$data
  D$headers <- TD$headers
  D$intype <- TD$type
  for(fname in colnames(D$input)){
    camel[[tolower(fname)]] <- fname
  }
  colnames(D$input) <- tolower(colnames(D$input))
  if('contig' %in% colnames(D$input)){
    D$input <- D$input %>% rename(seqname=contig)
    camel[['seqname']] <- 'Contig'
  }
  logWrite(paste('#LOAD',nrow(D$input),"input lines loaded from",filename))
  for(pfield in c("start","end","pos")){
    if(pfield %in% colnames(D$input)){
      D$input[[pfield]] <- as.integer(D$input[[pfield]])
    }
  }
}else{
  if(settings$rscript){
    quit("no",2)  
  }else{
    stop(paste("Cannot find input table:",filename))
  }
}

# seqlen=FILE : File of SeqName, SeqLen
# > D$seqlen : data frame of SeqName and SeqLen
#i# Needed for reverse complementing positions. Can be in the input table.
filename <- settings$seqlen
delimit <- "\t"
if(endsWith(filename,".txt")){
  delimit <- " "
}
if(endsWith(filename,".csv")){
  delimit <- ","
}
if(file.exists(filename)){
  D$seqlen <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="")
  logWrite(paste('#LOAD',nrow(D$seqlen),"seqlen lines loaded from",filename))
}else{
  if(filename != ""){
    stop(paste("Cannot find input table:",filename))
  }
}

# fields=DICT : list of field:fieldname, separated by commas (fieldname is in input)
# > settings$fields : list object
for(fname in names(settings$fields)){
  D$input <- rename(D$input,!!fname := settings$fields[[fname]])
}

    
#2# Transform input data using trandsb 
DF <- D
D <- DF$input
transdb <- DF$transdb

### ~ Optional splitting of data that is split by a transformation ~~~~~~~~~~~~~ ###
if(settings$split){
  #i# Establish which regions need splitting at beginning and/or end
  D$split5 <- 0
  D$split3 <- 0
  for(i in 1:nrow(transdb)){
    reg <- transdb[i,]
    if(reg$end < 1){
      if("seqlen" %in% names(DF) & reg$seqname %in% DF$seqlen$SeqName){
        reg$end <- DF$seqlen[DF$seqlen$SeqName == reg$seqname,]$SeqLen
      }else{
        reg$end <- maxlen
      }
    }
    #i# Construct a vector of TRUE or FALSE for D rows that need splitting
    D <- D %>%
      mutate(split5=ifelse(seqname == reg$seqname & start < reg$start & end > reg$start,reg$start,split5)) %>%
      mutate(split3=ifelse(seqname == reg$seqname & start < reg$end & end > reg$end,reg$end,split3))
  }
  nsplit <- sum(D$split5 > 0 | D$split3 > 0)
  logWrite(paste("#SPLIT",nsplit,"regions to be split at transformation boundaries (split=TRUE)"))
  #i# Construct split regions
  if(nsplit > 0){
    splitD <- 
      bind_rows(D %>% filter(split5 > 0) %>% mutate(end=split5-1),
        D %>% filter(split5 > 0 & split3 ==0) %>% mutate(start=split5),
        D %>% filter(split5 > 0 & split3 > 0) %>%  mutate(start=split5, end=split3),
        D %>% filter(split5 == 0 & split3 > 0) %>%  mutate(end=split3),
        D %>% filter(split3 > 0) %>%  mutate(start=split3+1))
    #i# Filter and combine original regions
    D <- bind_rows(D %>% filter(split5 == 0 & split3 == 0), splitD) %>%
      select(-split5,-split3)
  }
}

### ~ Join data ~~~~~~~~~~~~~~~~~~~~~~ ###
DF$output <- transformLocations(D,transdb,strict=settings$filter,logwrite=TRUE)
D <- DF

### ~ Output data ~~~~~~~~~~~~~~~~~~~~~~ ###
#!# Need to sort out CamelCase output and headers if required
fappend <- FALSE
if(length(D$headers) > 0){
  cat(paste(D$headers,collapse="\n"),file=settings$out)
  cat("\n",file=settings$out,append=TRUE)
  fappend <- TRUE
}

if(settings$camelcase){
  for(fname in names(camel)){
    if(fname %in% colnames(D$output)){
      #D$output <- rename(D$output,!!fname := camel[[fname]])
      D$output <- rename(D$output,!!camel[[fname]] := all_of(fname))
    }
  }
}

write.table(D$output,settings$out,sep="\t",quote=FALSE,row.names=FALSE,col.names=! D$intype %in% c("busco","gff"),append=fappend)
logWrite(paste('#SAVE',nrow(D$out),"transformed data output to",settings$out))


################# ::: CODE TO ADAPT ::: ######################


if(settings$rscript){
  quit("no",2)  
}

### ~ Reverse complement sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
revComp <- function(dnaseq){
    compvec <- c("A" = "1", "T" = "2", "G" = "3", "C" = "4", "a" = "5", "t" = "6", "g" = "7", "c" = "8", "u" = "6", "U" = "2",
                 "1" = "T", "2" = "A", "3" = "C", "4" = "G", "5" = "t", "6" = "a", "7" = "c", "8" = "g")
    newseq <- str_replace_all(dnaseq, compvec)
    newseq <- stringi::stri_reverse(newseq)
    return(newseq)
}

### ~ FASTA to sequence list (dictionary) ~~~~~~~~~~~~~~~~~~~~~~ ###
# > filename: fasta file
# > seqlist: list of name=sequence loaded from filename
# > shortname: whether to only use first word of name for seqlist name
fastaToList <- function(filename,seqlist=list(),shortname=TRUE){
  indata = readLines(filename)
  i = 1
  logWrite(paste0("Loading sequence data from ",filename,"..."))
  while(i < length(indata)){
    seqname = indata[i]
    if(startsWith(seqname,">")){
      cat(".", file = stderr())
      seqname = substr(seqname,2,nchar(seqname))
      if(shortname){
        seqname <- strsplit(seqname,split=" ",fixed=TRUE)[[1]][1]
      }
      seqlist[[seqname]] = indata[i+1]
      i = i + 2
    }else{
      seqlist[[names(seqlist)[length(seqlist)]]] <- paste0(seqlist[[names(seqlist)[length(seqlist)]]],seqname)
      i = i + 1
    }
  }
  cat("\n", file = stderr())
  logWrite(paste("#FASTA Sequence data for",length(seqlist),"sequences loaded from",filename))
  return(seqlist)
}

### ~ Sequence list (dictionary) to FASTA ~~~~~~~~~~~~~~~~~~~~~~ ###
# > filename: fasta file
# > seqlist: list of name=sequence loaded from filename
# > seqdesc: optional list of names to descriptions
listToFasta <- function(filename,seqlist,seqdesc=list(),append=FALSE){
  # Write output to filename
  for(seqname in sort(names(seqlist))){
    outname <- seqname
    if(seqname %in% seqdesc){
      outname <- paste(seqname, seqdesc[[seqname]])
    }
    cat(paste0(">",outname,"\n"), file=filename,append=append)
    append=TRUE
    cat(paste0(seqlist[[seqname]],"\n"), file=filename, append=TRUE)
  }
}


### ~ Remake sequence list from table ~~~~~~~~~~~~~~~~~~~~~~ ###
# > seqdb: tibble including NewSeqName, Chunk, SeqName, Start, End, Strand
# > seqlist: list of SeqName=sequence from which to build new sequences
# > makeD$seqlist: list of sequences to return
# > makeD$seqdb: updated tibble with NewStart and NewEnd
# > makeD$transdb: translation table matching the new sequences
makeSeqs <- function(seqdb,seqlist,newseqlist=list()){
  # Sort into the order of constructing new sequences
  newnames <- unique(seqdb[seqdb$SeqName %in% names(seqlist),]$NewSeqName)
  seqdb <- seqdb %>% filter(NewSeqName %in% newnames) %>% arrange(NewSeqName,Chunk)
  # End of sequence updates
  for(seqname in seqdb[seqdb$End < 1,]$SeqName){
    #print(seqname)
    seqlen <- nchar(seqlist[[seqname]])
    seqdb[seqdb$End < 1 & seqdb$SeqName == seqname,]$End <- seqlen
  }
  transdb <- seqdb
  seqdb$NewStart <- 0
  seqdb$NewEnd <- 0

  # Generate transdb
  transdb$Shift <- 0
  transdb$RevComp <- FALSE
  transdb[transdb$Strand == "-",]$RevComp <- TRUE
  # Make new sequences
  for(newseqname in newnames){
    ishift <- 0
    newseq <- ""
    chunks <- seqdb %>% filter(NewSeqName == newseqname)
    for(i in 1:nrow(chunks)){
        seqi <- seqdb$NewSeqName == newseqname & seqdb$Chunk == chunks$Chunk[i]
        seqdb[seqi,]$NewStart <- ishift + 1
        # - SeqName can be "Gap"
        if(chunks$SeqName[i] == "Gap"){
          newseq <- paste0(c(newseq,rep("N",chunks$End[i])),collapse="")
          ishift <- ishift + chunks$End[i]
        }else{
        # Otherwise, extract regions
          fullseq <- seqlist[[chunks$SeqName[i]]]
          partseq <- substr(fullseq,chunks$Start[i],chunks$End[i])
          if(chunks$Strand[i] == "-"){
            partseq <- revComp(partseq)
          }
          newseq <- paste0(newseq,partseq)
          transdb[seqi,]$Shift <- (ishift - chunks$Start[i] + 1)
          ishift <- ishift + chunks$End[i] - chunks$Start[i] + 1
        }
        seqdb[seqi,]$NewEnd <- ishift
        if(ishift != nchar(newseq)){
            writeLines(paste(newseqname,chunks$Chunk[i],"->",nchar(newseq),"bp but should be",ishift,"bp"))
        }
    }
    newseqlist[[newseqname]] <- newseq
  }
  return(list(seqlist=newseqlist,seqdb=seqdb,transdb=transdb))
}


### ~ Make transdb from newseq table ~~~~~~~~~~~~~~~~~~~~~~ ###
# > seqdb: tibble including NewSeqName, Chunk, SeqName, Start, End, Strand
# > seqlist: list of SeqName=sequence from which to build new sequences
# < transdb: Transformation table to return: SeqName, Start, End, NewSeqName, Shift, RevComp
makeTransDB <- function(seqdb,seqlist){
  # Sort into the order of constructing new sequences
  transdb <- seqdb %>% filter(SeqName %in% names(seqlist)) %>% arrange(desc(NewSeqName),Chunk)
  newnames <- unique(transdb$NewSeqName)
  transdb$Shift <- 0
  transdb$RevComp <- FALSE
  transdb[transdb$Strand == "-",]$RevComp <- TRUE
  # End of sequence updates
  for(seqname in transdb[transdb$End < 1,]$SeqName){
    #print(seqname)
    seqlen <- nchar(seqlist[[seqname]])
    transdb[transdb$End < 1 & transdb$SeqName == seqname,]$End <- seqlen
  }
  # Make new sequences
  for(newseqname in newnames){
    ishift <- 0
    chunks <- transdb %>% filter(NewSeqName == newseqname)
    for(i in 1:nrow(chunks)){
        # - SeqName can be "Gap"
        if(chunks$SeqName[i] == "Gap"){
          ishift <- ishift + chunks$End[i]
        }else{
        # Otherwise, extract regions
          treg <- transdb$NewSeqName == newseqname & transdb$Chunk == chunks$Chunk[i]
          transdb[treg,]$Shift <- (ishift - chunks$Start[i] + 1)
          ishift <- ishift + chunks$End[i] - chunks$Start[i] + 1
        }
    }
  }
  return(transdb)
}





