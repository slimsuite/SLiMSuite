########################################################
### GapMap.                                    ~~~~~ ###
### VERSION: 0.0.0                             ~~~~~ ###
### LAST EDIT: 27/05/22                        ~~~~~ ###
### AUTHORS: Richard Edwards 2022              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# Replace with one-line description of code purpose here.

####################################### ::: HISTORY ::: ############################################
# v0.0.0 : Initial version based on template.R
program <- "GapMap"
version <- "v0.0.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript template.R seqin=FILE depfile=FILE basefile=STR
# : pafin=FILE = PAF file to process [$BASEFILE.paf]
# : minmapq=INT = Minimum MapQ value to keep the alignment [60]
# : minctgpc=NUM = Minimum percentage of CtgLen to be covered by AlnLen [90]
# : busco=FILE = Full BUSCO results table for gap-filled assembly []
# : prevbusco=FILE = Full BUSCO results table for pre-gap-filled assembly []
# : basefile=STR = Prefix for file outputs [gapmap]
# : threads=NUM = Number of threads to use. [0 will use detectCores() - 1] (Not yet implemented)

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# os.popen('Rscript {0}deppurgehap.R {1}'.format(rdir, optionstr)).readlines()

####################################### ::: OUTPUTS ::: ############################################
#!# List the outputs here

####################################### ::: TODO ::: ############################################
# [ ] : Update the default inputs for BUSCO etc.
# [ ] : Add loading of gap table and filtering of remaining gaps from filled gaps.
#!# Add parallel processing of regions if speedup needed.

####################################### ::: SETUP ::: ############################################
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(pafin="",minmapq=60,minctgpc=90,
                busco="",prevbusco="",
                debug=FALSE,
                basefile="gapmap",rscript=TRUE,rdir="",
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
for(cmd in c("minmapq")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
for(cmd in c()){
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
if(settings$pafin %in% c("","None","none")){
  settings$pafin <- paste0(settings$basefile,".paf")
}
# LogWrite the PAFIN file
if(! file.exists(settings$pafin)){
  print("PAF file not found")
  #i# Print warning and quit
}

#i# Load in PAF file and add headers
pafdb <- read.table(settings$pafin,fill=TRUE,sep="\t",header=FALSE,row.names = NULL,quote="\"",comment.char="",stringsAsFactors=FALSE) %>% select(1:12)
colnames(pafdb) <- c("Ctg","CtgLen","CtgStart","CtgEnd","Strand","SeqName","SeqLen","Start","End","MatchNum","AlnLen","MapQ")
logWrite(paste('#PAF',nrow(pafdb),"PAF lines loaded from",settings$pafin))

#i# Need to filter by MapQ
pafdb <- pafdb %>% filter(MapQ >= settings$minmapq)
logWrite(paste('#PAF',nrow(pafdb),"PAF lines from",settings$pafin,"after filtering to MapQ >=",settings$minmapq))

pafdb$MatchPerc <- 100.0 * pafdb$MatchNum / pafdb$AlnLen
pafdb$CtgPerc <- 100.0 * pafdb$AlnLen / pafdb$CtgLen
pafdb <- pafdb %>% filter(CtgPerc >= settings$minctgpc)
logWrite(paste('#PAF',nrow(pafdb),"PAF lines from",settings$pafin,"after filtering to CtgPerc >=",settings$minctgpc))

# Split contig names into ScfName, ScfStart and ScfEnd
pattern <- "^(.*?)\\.(\\d+-\\d+)$"

# Extract the desired parts using str_extract()
#i# Not sure why this no longer produces ScfPos.1 and ScfPos.2
pafdb <- pafdb %>%
  mutate(
    ScfName = str_extract(Ctg, pattern,group=1),
    numbers = str_extract(Ctg, pattern,group=2),
    # Remove the period and split the numbers part
    ScfPos = str_remove(numbers, "\\.") %>% str_split(., "-", simplify = TRUE)
  ) %>% select(-numbers)
# Sort by ScfName and SeqName - gaps will be within a Scf (but want to identify spanning contigs later)
# Add clipping 5' and 3' of the alignment. This is 5' and 3' relative to SeqName not the mapped scaffold
pafdb$ScfStart <- pafdb$ScfPos[,1]
pafdb$ScfEnd <- pafdb$ScfPos[,2]
pafdb[is.na(pafdb$ScfName),]$ScfStart <- 1
pafdb[is.na(pafdb$ScfName),]$ScfEnd <- pafdb[is.na(pafdb$ScfName),]$CtgLen
pafdb <- pafdb %>% arrange(ScfName,SeqName,Start,End) %>%
  mutate(ClipStart=CtgStart-1,ClipEnd=CtgLen-CtgEnd)
pafdb[pafdb$Strand == "-",]$ClipStart <- pafdb[pafdb$Strand == "-",]$ClipEnd
pafdb[pafdb$Strand == "-",]$ClipEnd <- pafdb[pafdb$Strand == "-",]$CtgStart-1

#i# Split into scaffolds with ScfName for gap mapping and scaffolds with is.na(ScfName) for finding unplaced contigs mapped into gaps
scfdb <- pafdb %>% filter(! is.na(ScfName))
ctgdb <- pafdb %>% filter(is.na(ScfName))

# Add shifted SeqName and Start and then convert to gaps
gapdb <- scfdb %>% rename(Flank3=Start,Flank5=End)
gapdb$NextSeq <- gapdb$SeqName[c(2:nrow(gapdb),1)]
gapdb$Flank3 <- gapdb$Flank3[c(2:nrow(gapdb),1)]
gapdb$Clip5 <- gapdb$ClipEnd
gapdb$Clip3 <- gapdb$ClipStart[c(2:nrow(gapdb),1)]
gapdb <- gapdb %>% filter(SeqName==NextSeq) %>%
  mutate(EndDist = pmin(Flank5, SeqLen - Flank3), 
         EndPerc = 100.0 * EndDist / SeqLen, 
         Start=Flank5+1, End=Flank3-1, RegLen=Flank3-Flank5-1
  ) %>%
  mutate(
    Strand = if_else(Start < End, "+", "-"),
    Temp = if_else(Strand == "-", End, Start),
    End = if_else(Strand == "-", Start, End),
    Start = Temp
  ) %>%
  select(SeqName,Start,End,Strand,RegLen,SeqLen,EndDist,EndPerc,Clip5,Clip3) %>%
  mutate(FillLen=RegLen-Clip5-Clip3)

#!# NOTE: Remaining gaps need to be removed from this file, or flagged with a different tag.
# gabd <- gapdb %>% mutate(Type="gap")
#i# Load gap table
#!# regFeatures(regD,ftD,stranded=FALSE,overlap="any",flank=0,maxsize=100){

# Want to also add the start and end of each sequence
enddb <- pafdb %>%
  group_by(SeqName) %>%
  summarise(
    LowestStart = min(Start),
    HighestEnd = max(End),
    SeqLen = max(SeqLen)
  ) %>% mutate(Strand="+")

# Join this in
enddb <- bind_rows(# gapdb %>% mutate(Type="gap"),
  enddb %>% mutate(Start=1,End=LowestStart-1,Type="end5") %>% select(SeqName,Start,End,Strand,SeqLen,Type) %>% mutate(RegLen=End-Start+1),
  enddb %>% mutate(Start=HighestEnd+1,End=SeqLen,Type="end3") %>% select(SeqName,Start,End,Strand,SeqLen,Type) %>% mutate(RegLen=End-Start+1)
)

#!# gapdb and enddb can now be used to check with DepthKopy and Quenda.
outfile = paste(settings$basefile,"gap.tsv",sep=".",collapse=".")
logWrite(paste("#SAVE",nrow(gapdb),"Filled gap positions output to",outfile))
write.table(gapdb,outfile,sep="\t",quote=FALSE,row.names=FALSE)

outfile = paste(settings$basefile,"end.tsv",sep=".",collapse=".")
logWrite(paste("#SAVE",nrow(enddb),"End positions output to",outfile))
write.table(enddb,outfile,sep="\t",quote=FALSE,row.names=FALSE)



### ~~~~~~~~~~ Contig spanning and BUSCO analysis ~~~~~~~~~~~~~~~~~~~ ###
#!# Also want to look at ctgdb now and see if any of the contigs match to gap-filled regions
ctgspan <- full_join(gapdb, ctgdb %>% select(SeqName,Ctg,CtgLen,Start,End) %>% rename(CtgStartPos=Start,CtgEndPos=End),relationship = "many-to-many") %>% filter(End>=CtgStartPos,Start<=CtgEndPos)
logWrite(paste("#CTG",nrow(ctgspan),"contigs span filled gaps"))
if(nrow(ctgspan)>0){
  outfile = paste(settings$basefile,"ctgspan.tsv",sep=".",collapse=".")
  logWrite(paste("#SAVE",nrow(ctgspan),"contigs spanning gaps output to",outfile))
  write.table(ctgspan,outfile,sep="\t",quote=FALSE,row.names=FALSE)
}


if(file.exists(settings$busco)){
  busdb <- loadTable(settings$busco)$data
  buscospan <- full_join(gapdb %>% select(-Strand), busdb %>% select(Contig,BuscoID,Length,Start,End,Strand,Status,Score) %>% rename(SeqName=Contig,StartPos=Start,EndPos=End),relationship = "many-to-many") %>% filter(End>=StartPos,Start<=EndPos)
  logWrite(paste("#BUSCO",nrow(buscospan),"BuscoID span filled gaps"))

  if(file.exists(settings$prevbusco)){
    prevdb <- loadTable(settings$prevbusco)$data
    buscospan <- left_join(buscospan, 
                           prevdb %>% select(BuscoID,Status,Score,Contig,Start,End,Strand) %>% 
                             rename(PrevSeq=Contig,PrevStart=Start,PrevEnd=End,PrevStrand=Strand,PrevScore=Score,PrevStatus=Status)
    )
  }
  if(nrow(buscospan)>0){
    outfile = paste(settings$basefile,"buscospan.tsv",sep=".",collapse=".")
    logWrite(paste("#SAVE",nrow(buscospan),"BUSCOs spanning gaps output to",outfile))
    write.table(buscospan,outfile,sep="\t",quote=FALSE,row.names=FALSE)
  }
  
  #!# Do a quick comparison of the BUSCOs
  buscomp <- full_join(busdb %>% select(BuscoID,Status,Score) %>% group_by(BuscoID) %>% slice_max(Score) %>% unique(), 
                       prevdb %>% select(BuscoID,Status,Score) %>% group_by(BuscoID) %>% slice_max(Score) %>%
                         rename(PrevScore=Score,PrevStatus=Status) %>% unique()
  ) %>% filter(Status != PrevStatus) %>%
    mutate(Change=Score-PrevScore)
  
  
  #!# These changes are really interesting in the context of what HyPo is up to.
  #?# Are these gap-filling changes, or are they HyPo error correction
  summary(buscomp)
  logWrite(paste(nrow(buscomp),"BUSCO genes have changed classification"))
  table(buscomp$Status,buscomp$PrevStatus)
  saveTable(buscomp,suffix="busclass","BUSCO genes with altered classification")
  
  busdiff <- tibble(Status=c("Complete","Duplicated","Fragmented","Missing"))
  busdiff$Gapped <- 0
  busdiff$Filled <- 0
  for(i in 1:nrow(busdiff)){
    busdiff$Gapped[i] <- nrow(buscomp %>% filter(PrevStatus==busdiff$Status[i]))
    busdiff$Filled[i] <- nrow(buscomp %>% filter(Status==busdiff$Status[i]))
  }
  busdiff$GapFill <- busdiff$Filled - busdiff$Gapped
  saveTable(busdiff,suffix="busdiff","summary BUSCO Status shifts")
}



##### ======================== Tidy and Finish ======================== #####
if(file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE GapMap.R finished.")