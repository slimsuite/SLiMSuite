########################################################
### DepPurgeHap                               ~~~~~ ###
### VERSION: 0.1.0                             ~~~~~ ###
### LAST EDIT: 15/04/22                        ~~~~~ ###
### AUTHORS: Richard Edwards 2022              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for replacing Purge_haplotigs in Diploidocus

####################################### ::: HISTORY ::: ############################################
# v0.0.0 : Initial version based on depthcopy.R.
# v0.1.0 : Initial working version.
program <- "DepPurgeHap"
version <- "v0.1.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript deppurgehap.R depfile=FILE [busco=FILE] [scdepth=INT] [gablam=PREFIX] 
# : depfile=FILE = Full depth data table (one value per base, pseudo-fasta format)
# : busco=FILE = Full BUSCO table. If a number, will interpret as scdepth value.
# : gablam=PREFIX = Prefix for GABLAM output (*.gablam.tdt and *.Local.tdt)
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


### ~ Load Depth File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# seqVector(filename,intvec=TRUE) - Returns list of seqname=vector
#i# loadTable(filename) 

#i# To filter to subset of sequences:
#i# filterTableSeq(D,seqnames,seqfield="SeqName",logid="#FILTER") - Returns filtered tibble/dataframe

### ~ Load BUSCO File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load BUSCO table into tibble
#i# buscodb = buscoTable(filename)
buscoTable <- function(filename,status="Complete"){
  TD <- loadTable(filename) 
  buscodb <- TD$data
  buscodb = buscodb[buscodb$Status %in% status,]
  logWrite(paste('#BUSCO',nrow(buscodb),paste0(status,collapse="+"),"BUSCO genes loaded from",filename))
  return(buscodb)
}
#i# dupdb = buscoDupTable(filename)
buscoDupTable <- function(filename){
  buscodb = buscoTable(filename,"Duplicated")
  return(buscodb)
}


### ~ Calculate SC Depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Straight mode of depth vector
getmode <- function(v) {
  uniqv <- unique(v)
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}
#i# Density plotting function
densModePlot <- function(bdens,centre,scdepth,plotmain){
  
  plot(bdens,xlim=c(0,scdepth*2), lwd = 2,col="blue",
       main=plotmain,ylab="Depth Density",xlab="XDepth",
       ylim=c(0,max(bdens$y)))
  
  abline(v=centre,col="grey")
  text(centre-2,0,round(centre,2),adj=c(1,0),col="grey")
  
  dmode=bdens$x[bdens$y == max(bdens$y)]
  abline(v=dmode,col="steelblue")
  text(dmode+2,0,round(dmode,2),adj=0,col="steelblue")
  abline(h=max(bdens$y),col="steelblue")
  
  abline(v=scdepth,col="red")
  text(scdepth+2,max(bdens$y)/2,round(scdepth,2),adj=0,col="red")
  
}
#i# Density-based mode calculation
densModeZoom <- function(depvec,adjust=16,plotbase=NA,plotmain=NA){
  depvec = depvec[!is.na(depvec)]
  if(length(depvec) < 1){
    if(settings$debug){
      logWrite("DepVec length <1")
    }
    return(0)
  }
  centre = getmode(depvec)
  if(centre <= 1){
    if(settings$debug){
      logWrite("Centre <= 0")
    }
    return(0)
  }
  #newvec = depvec[depvec <= (centre*4) & depvec >= centre/4]
  newvec = depvec[depvec <= max(1000,centre*4)]
  if(length(newvec) < 1){
    if(settings$debug){
      logWrite("NewVec length <1")
    }
    return(0)
  }
  #!# Add a temporary fudge for vary sparse and high copy regions
  mindensn <- 10
  if(length(newvec) < mindensn){
    if(settings$debug){
      logWrite('#DENS Fudging density calculation')
    }
    return(mean(depvec))
  }
  n = 2048
  while(max(newvec)*5 > n){
    n = n * 2
  }
  depdens = density(newvec,n=n,adjust=adjust)
  scdepth = depdens$x[depdens$y == max(depdens$y)]
  return(scdepth)
}

### ~ Predict CN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Extract depth profiles from BUSCO Complete data
buscoDepData <- function(buscofile,deplist){
  buscodb = buscoTable(buscofile)
  buscovec = list()  # List of BuscoID: depvec
  depvec = c()       # Total depth vector
  # Extract depth vectors
  for(i in 1:nrow(buscodb)){
    seqname = buscodb$Contig[i]
    posi = buscodb$Start[i]
    posj = buscodb$End[i]
    seqdep = deplist[[seqname]][posi:posj]
    buscovec[[buscodb$BuscoID[i]]] = seqdep
    if(settings$debug){
      writeLines(paste(i,buscodb$BuscoID[i],seqname,posi,posj,length(seqdep)))
    }
    depvec = c(depvec,seqdep)
    cat(".", file = stderr())
  }
  cat("\n", file = stderr())
  writeLines("Calculating BUSCO SC Depth...")
  print(summary(depvec))
  #i# Adjust when excess zero-coverage regions are messing things up
  depvec = depvec[!is.na(depvec)]
  centre = getmode(depvec)
  if(centre < 1){
    logWrite("#WARN Warning: pure BUSCO modal depth <1X - trimming bases <1X coverage")
    depvec = depvec[depvec >=1 ]
  }
  #i# Calculate SCDepth from BUSCO genes
  scdepth = densModeZoom(depvec,adjust=settings$adjust,plotbase=settings$basefile,plotmain="BUSCO")
  #!# Generate PNG BUSCO depth profile
  logWrite(paste("#SCDPETH BUSCO SC Depth:",round(scdepth,2)))
  ### ~ Save SC Depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  if(scdepth > 0){
    outfile = paste(settings$depfile,"scdepth",sep=".",collapse=".")
    logWrite(paste(nrow(buscodb),"BUSCO Complete SC Depth output to",outfile))
    write(scdepth,outfile)
    if(settings$scdepth > 0){
     if(settings$buscocn){
        logWrite(paste("#SCDPETH Over-riding scdepth=NUM SC Depth:",round(scdepth,2)))
     }else{
      scdepth = settings$scdepth
      logWrite(paste("#SCDPETH Using scdepth=NUM SC Depth:",round(scdepth,2)))
     }
    }
  }else{
    logWrite(paste(nrow(buscodb),"BUSCO Complete SC Depth calculation failed!"))
    if(settings$scdepth > 0){
      scdepth = settings$scdepth
      logWrite(paste("#SCDPETH Using scdepth=NUM SC Depth:",round(scdepth,2)))    
    }else{
      quit("no",1)
    }
  }
  return( list(depvec=depvec, buscodb=buscodb, buscovec=buscovec, scdepth=scdepth) )
}

#i# Combine BUSCO data with depth profiles
#># buscodat = buscoDepCalc(buscodat)
buscoDepCalc <- function(buscodat){
  depvec = buscodat$depvec
  scdepth = buscodat$scdepth
  buscodb = buscodat$buscodb
  buscovec = buscodat$buscovec
  if(nrow(buscodb)<1){
    logWrite('Zero BUSCO rows: skipping Depth Calculation')
    return( list(depvec=depvec, buscodb=buscodb, buscovec=buscovec, scdepth=scdepth) )
  }
  # Setup new fields. Some of these are just for testing
  buscodb$MeanX = 0
  buscodb$MedX = 0  
  buscodb$ModeX = 0  
  buscodb$DensX = 0  
  buscodb$CN = 0
  buscodb$CIsyst = 0.0
  buscodb$CIrand = 0.0
  # Calculate gene stats
  for(i in 1:nrow(buscodb)){
    #logWrite(buscodb$BuscoID[i])
    seqdep = buscovec[[buscodb$BuscoID[i]]]
    if(length(seqdep) >= 1){
      buscodb$MeanX[i] = mean(seqdep)
      buscodb$ModeX[i] = getmode(seqdep)
      buscodb$MedX[i] = median(seqdep)
      buscodb$DensX[i] = densModeZoom(seqdep,adjust=settings$adjust)
    }else{
      logWrite(paste("#WARN BUSCOID",buscodb$BuscoID[i],"has zero-length depth vector"))
    }
    cat(".", file = stderr())
  }
  cat("\n", file = stderr())
  # Return updated data
  logWrite("#BUSCO BUSCO depth calculations complete.")
  # Add Homology / Kmer data if present
  buscodb$DensK = NA
  buscodb$SelfK = NA
  buscodb$HomPC = NA
  writeLines("Calculating BUSCO kmers/homology...")
  # Calculate region stats
  prevname = ""
  for(i in 1:nrow(buscodb)){
    seqname = buscodb$Contig[i]
    if(seqname != prevname){
      prevname = seqname
      if(settings$debug){
        cat(seqname)
      }
    }
    cat(".", file = stderr())
    seqi = buscodb$Start[i]
    seqj = buscodb$End[i]
    if(file.exists(katfile)){
      #i# katlist = raw read kmer frequencies
      seqdep = depVector(seqname,seqi,seqj,katlist)
      buscodb$DensK[i] = densModeZoom(seqdep,adjust=settings$adjust)
    }
    if(file.exists(katself)){
      #i# katself = assembly kmer frequencies
      seqdep = depVector(seqname,seqi,seqj,selflist)
      buscodb$SelfK[i] = densModeZoom(seqdep,adjust=settings$adjust)
    }
    if(file.exists(homfile)){
      #i# homlist = Assembly vs self homology
      seqdep = depVector(seqname,seqi,seqj,homlist)
      #buscodb$HomX[i] = densModeZoom(seqdep,adjust=settings$adjust)
      buscodb$HomPC[i] = 100.0 * sum(seqdep > 0) / length(seqdep)
    }
  }
  cat("\n", file = stderr())
  logWrite("#BUSCO BUSCO kmer/homology calculations complete.")
  return( list(depvec=depvec, buscodb=buscodb, buscovec=buscovec, scdepth=scdepth))
}


################################## ::: RUN CODE ::: #######################################

##### ======================== Report key inputs ======================== #####

### ~ Key inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Depth file of read depths per base
depfile = settings$depfile
logWrite(paste("Depth File:",depfile))
if(! file.exists(depfile)){
  quit("no",2)  
}
#i# BUSCO table or scdepth
buscofile = settings$busco
scdepth = as.numeric(settings$scdepth)
if(file.exists(buscofile)){
  logWrite(paste("BUSCO Full File:",buscofile))
}else{
  logWrite(paste("SC Depth:",round(scdepth,2)))
}
logWrite('#RCODE Setup complete.')


##### ======================== Load data ======================== #####

### ~ Load Depth data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# 1. Load the depth data list
deplist = seqVector(depfile)
settings$seqnames = names(deplist)


### ~ Load GABLAM tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
locfile <- paste0(settings$gablam,".local.tdt")
gabfile <- paste0(settings$gablam,".gablam.tdt")
gabdb <- tibble()
locdb <- tibble()
if(file.exists(locfile) & file.exists(gabfile)){
  gabdb <- loadTable(gabfile,delimit="ext")
  gabdb <- gabdb$data
  locdb <- loadTable(locfile,delimit="ext")
  locdb <- locdb$data
}else{
  logWrite(paste0('#GABLAM ',settings$gablam,".* GABLAM files not found."))
}


### ~ Calculate SC Depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# 2. Calculate SC read depth if mot provided
buscoMean = 0
buscoSD = 0
scdepth <- settings$scdepth
if(scdepth < 1 & file.exists(buscofile)){
  #!# Add ability to re-use loaded BUSCO data and given scdepth=NUM value
  buscodat = buscoDepData(buscofile,deplist)
  scdepth = buscodat$scdepth
}

### ~ Calculate Sequence Coverage Bins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
phlow <- as.integer(scdepth / 4)
dupdep <- scdepth / 2
phmid <- as.integer(1.5 * dupdep)
phhigh <- as.integer(scdepth * 2 + 0.5) 
logWrite(paste0('#PURGE LowDep=',phlow,"X; MidDep=",phmid,"X; HighDep=",phhigh,"X"))
#i# 'SeqName','SeqLen','LowPerc','HapPerc','DipPerc','HighPerc','TopHit','SecHit','TopHitCov','MaxHitCov','PurgeHap','TopNum','SecNum'
#i# Generate initial table from depth list, then add details from the local and gablam table
hapdb <- tibble(SeqName=names(deplist))
for(field in c('SeqLen','LowPerc','HapPerc','DipPerc','HighPerc')){ #,'TopHit','SecHit','TopHitCov','MaxHitCov','PurgeHap','TopNum','SecNum')){
  hapdb[[field]] <- 0  
}
#i# Percentage depth bins
writeLines("Calculating depth bins...")
for(i in 1:nrow(hapdb)){
  cat(".", file = stderr())
  seqname <- hapdb$SeqName[i]
  seqdep <- deplist[[seqname]]
  hapdb$SeqLen[i] <- length(seqdep)
  hapdb$LowPerc[i] <- 100 * sum(seqdep < phlow) / hapdb$SeqLen[i]
  hapdb$HapPerc[i] <- 100 * sum(seqdep >= phlow & seqdep < phmid) / hapdb$SeqLen[i]
  hapdb$DipPerc[i] <- 100 * sum(seqdep >= phmid & seqdep <= phhigh) / hapdb$SeqLen[i]
  hapdb$HighPerc[i] <- 100 * sum(seqdep > phhigh) / hapdb$SeqLen[i]
}
cat("\n", file = stderr())
logWrite("#DEPBIN Depth bin calculations complete.")

#i# Top Hit and Second Hit
#i# Filter self-hits
gabdb <- gabdb[gabdb$Qry != gabdb$Hit,]
locdb <- locdb[locdb$Qry != locdb$Hit,]

backdb <- locdb

locdb <- backdb %>% arrange(Qry,Hit,QryStart,QryEnd) %>% select(-cs)
locdb$Length <- locdb$QryEnd - locdb$QryStart + 1
locdb$Pid <- locdb$Identity / locdb$Length
locdb$Keep <- TRUE

for(i in 2:nrow(locdb)){
  loci <- locdb[i-1,]
  locj <- locdb[i,]
  if(locj$QryStart > loci$QryEnd | loci$Qry != locj$Qry | loci$Hit != locj$Hit){ next }
  cat(".")
  #i# Combine overlap
  overlap <- (loci$QryEnd - locj$QryStart + 1) / 2
  if(locj$QryEnd > loci$QryEnd){
    idbp <- ((loci$Length - overlap) * loci$Pid) + ((locj$Length - overlap) * locj$Pid)
  }else{
    idbp <- ((loci$Length - overlap) * loci$Pid) + (overlap * locj$Pid)
    locdb$QryEnd[i] <- loci$QryEnd
  }
  locdb$Keep[i-1] = FALSE
  locdb$QryStart[i] <- loci$QryStart
  locdb$Length[i] <- locdb$QryEnd[i] - locdb$QryStart[i] + 1
  locdb$Identity[i] <- idbp
  locdb$Pid[i] <- (idbp / locdb$Length[i])
  if(locdb$Pid[i] > 1){
    print(locdb[(i-1):i,])
  }
}
cat("\n")
summary(locdb)

logWrite(paste("#LOCAL",nrow(locdb),"local hits merged to",sum(locdb$Keep)))
locdb <- locdb[locdb$Keep,]
locdb$HitCov <- 100 * (locdb$Length / locdb$QryLen) * locdb$Pid

#i# Compress per hit
matchdb <- locdb %>% rename(SeqName=Qry) %>% group_by(SeqName,Hit) %>% 
  summarise(SeqName=first(SeqName), Hit=first(Hit), TopHitCov=sum(HitCov))
#i# Filter on minhit
hitn <- nrow(matchdb)
matchdb <- matchdb[matchdb$TopHitCov >= settings$minhit,]
logWrite(paste("#HITCOV",hitn,"Qry-Hit pairs filter to",nrow(matchdb),"on minhit >=",settings$minhit))

#i# Compress to the best hits
maxdb <- matchdb %>% group_by(SeqName) %>% summarise(SeqName=first(SeqName),MaxHitCov=sum(TopHitCov))
topdb <- matchdb %>% group_by(SeqName) %>% slice_max(TopHitCov,n=1) %>% slice_head(n=1)
secdb <- matchdb %>% group_by(SeqName) %>% slice_max(TopHitCov,n=2) %>% slice_head(n=2) %>% slice_min(TopHitCov) %>% slice_tail(n=1) 

#i# Now combine the relevant summary stats with hapdb using left_join
purgedb <- left_join(hapdb,topdb %>% rename(TopHit=Hit))
purgedb <- left_join(purgedb,secdb %>% rename(SecHit=Hit) %>% select(SeqName,SecHit))
purgedb <- left_join(purgedb,maxdb)
purgedb <- left_join(purgedb,topdb %>% rename(Hap=SeqName,SeqName=Hit) %>% group_by(SeqName) %>% summarise(SeqName=first(SeqName),TopNum=n()))
purgedb <- left_join(purgedb,secdb %>% rename(Hap=SeqName,SeqName=Hit) %>% group_by(SeqName) %>% summarise(SeqName=first(SeqName),SecNum=n()))

#!# Fill in the missing data
purgedb[is.na(purgedb$TopHit),]$TopHit <- ""
purgedb[is.na(purgedb$SecHit),]$SecHit <- ""
purgedb[is.na(purgedb$TopHitCov),]$TopHitCov <- 0
purgedb[is.na(purgedb$MaxHitCov),]$MaxHitCov <- 0
purgedb[is.na(purgedb$TopNum),]$TopNum <- 0
purgedb[is.na(purgedb$SecNum),]$SecNum <- 0

#!# Add Rating
purgedb$PurgeHap <- "GABLAM"

#i# Output data
outfile <- paste0(settings$basefile,".deppurgehap.tsv")
write.table(purgedb, outfile, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
logWrite(paste("#SAVE DepPurgeHap data output to",outfile))
logWrite("#NOTE DepPurgeHap may have reciprocal TopHits for shorter contigs and repeats - additional filtering may be required.")

# 
# # ==> tiger.wtdbg2v1.racon2.10x.pilon2.scaffolds.purge.reassignments.tsv <==
# # #reassigned_contig      top_hit_contig  second_hit_contig       best_match_coverage     max_match_coverage      reassignment
# # Scaff10x_1000   Scaff10x_562    Scaff10x_676    98.37   918.99  REPEAT
# purgedb = db.addTable(purge,mainkeys=['#reassigned_contig'],expect=True,name='purge')
# purgedb.renameField('#reassigned_contig','SeqName')
# purgedb.renameField('top_hit_contig','TopHit')
# purgedb.renameField('second_hit_contig','SecHit')
# purgedb.renameField('best_match_coverage','TopHitCov')
# purgedb.renameField('max_match_coverage','MaxHitCov')
# purgedb.renameField('reassignment','PurgeHap')
# purgedb.index('TopHit')
# purgedb.index('SecHit')
# purgedb.addFields(['TopNum','SecNum'],evalue=0)
# for entry in purgedb.entries():
#   if entry['SeqName'] in purgedb.index('TopHit'): entry['TopNum'] = len(purgedb.index('TopHit')[entry['SeqName']])
# if entry['SeqName'] in purgedb.index('SecHit'): entry['SecNum'] = len(purgedb.index('SecHit')[entry['SeqName']])
# joinlist.append((purgedb,'SeqName'))



### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE DepPurgeHap.R finished.")
if(settings$rscript){
  quit("no",0)
}
