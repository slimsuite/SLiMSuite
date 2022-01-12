########################################################
### RJE Genomics Plot Functions                ~~~~~ ###
### VERSION: 0.3.0                             ~~~~~ ###
### LAST EDIT: 03/01/21                        ~~~~~ ###
### AUTHORS: Richard Edwards 2016              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for genomic plotting functions to be used by other scripts, 
# e.g. PAGSAT.R and samtools.R

################# ::: HISTORY ::: ######################
# v0.0.0 : Initial version based on PAGSAT v0.7.4 for use with samtools.R v0.0.0.
# v0.1.0 : Added SAMPhaser methods for used with SAMPhaser v0.4.0.
# v0.2.0 : Added Genome feature setup and plotting methods.
# v0.3.0 : Added some assembly functions.


################# ::: TO DO ::: ######################
# [Y] : Add function for transforming genomic locations (newseqname,shift,revcomp).
# [Y] : Add functions for generating new assembly from old assembly and SynBad-style edit table.


################# ::: SETUP ::: ######################

if(! "settings" %in% ls()){
  settings <- list(outlog=stdout())
}

if(! "logWrite" %in% ls()){
    if(! outlog %in% names(settings)){ settings$outlog <- stdout() }
    logWrite <- function(logstr){
      writeLines(paste0("[",date(),"] ",logstr),con=settings$outlog)
    }
}

################# ::: ASSEMBLY FUNCTIONS ::: ######################

### ~ Transform genomic locations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# > D : data frame with SeqName, Start, End, and Strand
# > transdb : data frame with SeqName, Start, End, NewSeqName, Shift, RevComp
# > strict=FALSE : boolean as to whether to exclude any regions in D that are not in transdb
# > logwrite=TRUE : Whether to summarise to logWrite()
# [ ] : Future > split=TRUE : whether to split partially overlapping regions rather than dumping them
transformLocations <- function(D,transdb,strict=FALSE,logwrite=TRUE){
    # Check that unique(SeqName, Start, End) is unique
    transdb <- transdb %>% filter(! SeqName == "Gap")
    ucheck <- transdb %>% select(SeqName, Start, End)
    ucheck <- unique(ucheck)
    if(nrow(ucheck) != nrow(transdb)){
       abort("transformLocations() was given transdb with non-unique Seqname, Start, End data")
    }
    maxlen <- max(D$End)
    if(max(transdb$End) < 1){
      transdb[transdb$End < 1,]$End <- maxlen
    }
    # Work through D and shift
    newD <- D
    newD$Edit <- FALSE
    for(i in 1:nrow(transdb)){
        reg <- transdb[i,]
        if(reg$SeqName == reg$NewSeqName & reg$Shift == 0 & ! reg$RevComp){ next }
        newDreg <- D$SeqName == reg$SeqName & D$Start >= reg$Start & D$End <= reg$End
        if(sum(newDreg) == 0){ next }
        if(sum(newD[newDreg,]$Edit) > 0){
            abort("transformLocations() is trying to edit same region multiple times")
        }
        writeLines(paste(sum(newDreg),"transformations"))
        newD[newDreg,]$SeqName <- reg$NewSeqName
        if(reg$RevComp){
            # Need to move everything to be from reg$End not 1 and swap
            newD[newDreg,]$Start <- reg$End - D[newDreg,]$End + 1 + reg$Shift
            newD[newDreg,]$End <- reg$End - D[newDreg,]$Start + 1 + reg$Shift
            if(sum(newDreg & D$Strand == "+")){ #   ("+" %in% newD[newDreg,]$Strand){
              newD[newDreg & D$Strand == "+",]$Strand <- "-"
            }
            if(sum(newDreg & D$Strand == "-")){ #if("-" %in% newD[newDreg,]$Strand){
              newD[newDreg & D$Strand == "-",]$Strand <- "+"
            }
        }else{
            # Otherwise, simple shift
            newD[newDreg,]$Start <- D[newDreg,]$Start + reg$Shift
            newD[newDreg,]$End <- D[newDreg,]$End + reg$Shift
        }
        newD[newDreg,]$Edit <- TRUE
    }
    # Summarise
    if(logwrite & "logWrite" %in% ls()){
      logWrite(paste(sum(newD$Edit),"of",nrow(D),"regions mapped to new names/locations"))
    }
    return(newD)
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


################# ::: PLOTTING FUNCTIONS ::: ######################

### ~ Generic Chromosome Plot Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# > xmin=0 : Min chromosome position WITHOUT scaling
# > xmax=1e5 : Max chromosome position WITHOUT scaling
# > ymin=0,ymax=1.0,ylab="",main="Chromosome Plot",ygrid=c(0.1,0.5)
chromPlotSetup = function(xmin=1,xmax=1e5,ymin=0,ymax=1.0,ylab="",main="Chromosome Plot",ygrid=c(0.1,0.5),yaxt="s"){
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  xmin = xmin/scaling
  xmax = xmax/scaling
  xlab=paste("Chromosome Position",units)
  plot(c(xmin,xmax),c(ymin,ymax),col="red",type="n",ylab=ylab,main=main,xlab=xlab,xaxs="i",yaxt=yaxt)
  plotGridlines(0,xmax,ymin,ymax,10000/scaling,50000/scaling,ygrid[1],ygrid[2])  
}


# SAMTools PLOTS

### ~ Chromosome Xcoverage plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# This plots depth of coverage across a chromosome or contig
# > chromdata = data table that has:
#   - $Chrom = Chromosome/Contig Locus Identifer
#   - $Pos = Position along chromosome
#   - #X = XDepth of coverage at position $Pos
# > chrom = Locus from chromdata to plot
# > onex = the X coverage for single copy (e.g. median X coverage) 
chromDepthPlot = function(chromdata,chrom,onex,dirplot=FALSE,dirdata=data.frame(),partplot=FALSE,ylab="XCoverage"){
  # Setup plot region and scaling parameters 
  scaling = settings$scaling
  plotdata = chromdata[chromdata$Locus==chrom,]
  plotdata = plotdata[(order(plotdata$Pos)),]
  ymax = max(plotdata$X,onex)
  ymax = min(ymax,10*onex)
  # Setup Empty Plot Region and Labels. 
  main = paste(chrom,"Read coverage")
  chromPlotSetup(xmax=max(plotdata$Pos),ymax=ymax,ylab=ylab,main=main,ygrid=c(onex/8,onex/2))
  
  # Plot Data
  if(dirplot==TRUE){
    dirplot = dirdata[dirdata$Locus==chrom,]
    dirplot = dirplot[(order(dirplot$Pos)),]
    lines(dirplot$Pos/scaling,dirplot$Len5,col="purple",type="l",lwd=1.5)
    lines(dirplot$Pos/scaling,dirplot$Len3,col="green",type="l",lwd=1.5)
  }
  if(partplot==TRUE){
    lines(plotdata$Pos/scaling,plotdata$FullX,col="lightblue",type="l",lwd=2)
    lines(plotdata$Pos/scaling,plotdata$PartX,col="red",type="l",lwd=2)
  }
  lines(plotdata$Pos/scaling,plotdata$X,col="blue",type="l",lwd=2)
  abline(h=onex,col="red")
}

#i# Function to add coverage of specific regions to a depth plot
#!# Assumes chromDepthPlot() already called with relevant Xdepth data
# > checkdata = data table that has:
#   - Locus	Start	End	
#   - Optional: Desc	
#   - Various SpanX, matching	spans (below) with at least $Span0
#   - MeanX = mean X depth across region
# > chrom = Chromosome/Locus for output
# > spans = vector of flank lengths for spanning read coverage
# > spancol = span colour - region colour will be transparent version
#!# Change to two colours for no/some spanning coverage, default red and blue: rgb(1,0,0,0.3) & rgb(0,0,1,0.3)
regionDepth = function(checkdata,chrom,spans=c(0),col=c(rgb(1,0,0),rgb(0,0,1))){
  scaling = settings$scaling
  chromdata = checkdata[checkdata$Locus==chrom,]
  if(dim(chromdata)[1] < 1){ return() }
  for(i in 1:dim(chromdata)[1]){  # Each row of check data
    plotdata = chromdata[i,]
    spancol = col[1]
    if(plotdata$Span0 > 0){
      spancol = col[2]
    }
    regcol = paste(spancol,"30",sep="")
    rect(plotdata$Start/scaling,-1,plotdata$End/scaling,plotdata$MeanX,lwd=2,col=regcol)  # Add transparent colour?
    for(flank in spans){
      field = paste("Span",flank,sep="")
      x = (plotdata$Start - flank) / scaling
      y = (plotdata$End + flank) / scaling
      h = plotdata[[field]]
      lines(c(x,y),c(h,h),lwd=2,col=spancol)
    }
  }

}





# PAGSAT PLOTS - Add more comments and descriptions 
#i# These can be found in pagsat-methods.R



# SAMPhaser PLOTS

################# ::: SAMPhaser PLOTTING FUNCTIONS ::: ######################

### ~ Plot a region and its xcoverage on an existing plot space ~~~~~~~~~~~~~~~~~ ###
#i# Plot area set up by samphaserStackPlot()
#i# hapdata: Locus	Block	Track	Start	End	HapStart	HapEnd	Haplotig	nSNP	nRID	MeanX

#># trackdepdata: Locus, Pos, X - where Locus is the Haplotype AccNum
#># trackhapdata: Locus	Block	Track	Start	End	HapStart	HapEnd	Haplotig HapAcc	nSNP	nRID	MeanX
#># onex: Median X depth
#># xshift = length to add to Start/End positions
#># ymin and ymax determine the vertical plot positions so that this code can be reused in different plots
#># tcol: colour for track
#i# The region is plotted at ymin+0.25onex to ymin+0.5onex; Main plot zero is at 0.5onex
regionXPlot = function(trackdepdata,trackhapdata,onex,ymin,ymax,tcol){
  scaling = settings$scaling
  # Set up local data to plot
  (xshift = trackhapdata$Start[1] - 1)
  # Plot the track
  rect(trackhapdata$Start/scaling,ymin+0.25*onex,trackhapdata$End/scaling,ymin+0.5*onex,col=tcol) 
  if(trackhapdata$HapStart>0){
    points(trackhapdata$HapStart/scaling,ymin+onex,pch=23,bg=tcol)
    points(trackhapdata$HapEnd/scaling,ymin+onex,pch=23,bg=tcol)
  }
  # Plot the Xcov data
  lines((trackdepdata$Pos+xshift)/scaling,trackdepdata$X+ymin+(0.5*onex),col="blue",type="l",lwd=2)
  abline(h=ymin+(1.5*onex),col="red")
}

regionXPlot2 = function(trackdepdata,trackhapdata,onex,ymin,tcol){
  scaling = settings$scaling
  # Set up local data to plot
  (xshift = trackhapdata$Start[1] - 1)
  # Plot the track
  rect(trackhapdata$Start/scaling,ymin,trackhapdata$End/scaling,ymin+onex/8,col=paste(tcol,"95",sep="")) 
  if(trackhapdata$HapStart>0){
    points(trackhapdata$HapStart/scaling,onex,pch=23,bg=paste(tcol,"95",sep=""))
    points(trackhapdata$HapEnd/scaling,onex,pch=23,bg=paste(tcol,"95",sep=""))
  }
  # Plot the Xcov data
  lines((trackdepdata$Pos+xshift)/scaling,trackdepdata$X,col=tcol,type="l",lwd=3)
}

# This plots a set of contigs or chromosomes against a query using local data and ctype depth data
#># depdata: Locus, Pos, X
#># hapdepdata: Locus, Pos, X - where Locus is the Haplotype AccNum
#># hapdata: Locus	Block	Track	Start	End	HapStart	HapEnd	Haplotig HapAcc	nSNP	nRID	MeanX
#># locus: Locus for plot
#># onex: Median X depth
#># ylim: Maximum multiplier of onex for ymax. (Avoids crazy depths messing things up)
samphaserStackPlot = function(depdata,hapdepdata,hapdata,locus,onex,ylim=10){
  ## Setup Empty Plot Region and Labels. chromPlotSetup will scale x-axis automatically.
  #i# The individual depth plot uses ygrid=c(onex/8,onex/2) and ymax=min(ymax,10*onex)
  #i# The region itself should be plotted at -0.25 to 0, with a gap of -0.5 to -0.25 between each A/B/C section
  (xcovmax = max(depdata[depdata$Locus==locus,]$X))
  (xmax = max(depdata[depdata$Locus==locus,]$Pos))
  # First establish onex value that covers full depth
  ymax = 2 * onex
  while(ymax < xcovmax){ ymax = ymax + onex/2 }
  # Add extra for the plot and gap
  (ymax = ymax + onex/2)
  ymax = min(ymax,ylim*onex)
  main = paste(locus,"haplotig read depth")
  chromPlotSetup(xmax=xmax,ymax=ymax*3,ylab="",main=main,ygrid=c(onex/8,onex/2),yaxt="n")
  
  ## Generate haplotig depth plots
  ymin = 0
  trackcol = list(A=soton$col[1],B=soton$col[2],C=soton$col[3])
  for(track in c("C","B","A")){
    text(0,(ymax+ymin+ymin)/2,track,pos=2,xpd=TRUE,cex=0.8,col=trackcol[[track]])
    for(hapacc in hapdata[hapdata$Locus==locus & hapdata$Track==track,]$HapAcc){
      (trackdepdata = hapdepdata[hapdepdata$Locus==hapacc,])
      (trackhapdata = hapdata[hapdata$HapAcc==hapacc,])
      regionXPlot(trackdepdata,trackhapdata,onex,ymin,ymin+ymax,trackcol[[track]])
    }
    ymin = ymin + ymax
  }
}

samphaserStackPlot2 = function(depdata,hapdepdata,hapdata,hapsnpdata,locus,onex,zoom=c(0,0),ylim=10){
  ## Setup Empty Plot Region and Labels. chromPlotSetup will scale x-axis automatically.
  #i# The individual depth plot uses ygrid=c(onex/8,onex/2) and ymax=min(ymax,10*onex)
  #i# The region itself should be plotted at -0.25 to 0, with a gap of -0.5 to -0.25 between each A/B/C section
  scaling = settings$scaling
  plotdata = depdata[depdata$Locus==locus,]
  (xcovmax = max(depdata[depdata$Locus==locus,]$X))
  (xmax = max(depdata[depdata$Locus==locus,]$Pos))
  xmin = 1
  #i# Zoom is already scaled. To do: Enable negative numbers to zoom from 3' end.
  if(zoom[1] > 0){ xmin = scaling * zoom[1]}
  if(zoom[2] > 0){ xmax = scaling * zoom[2]}
  # First establish onex value that covers full depth
  ymax = onex
  while(ymax < xcovmax){ ymax = ymax + onex/2 }
  ymax = min(ymax,ylim*onex)
  # Add extra for the plot and gap
  main = paste(locus,"haplotig read depth")
  chromPlotSetup(xmin=xmin,xmax=xmax,ymin=-onex/2,ymax=ymax,ylab="",main=main,ygrid=c(onex/8,onex/2))
  
  ## Generate haplotig depth plots
  ymin = -onex/2
  trackcol = list(A=soton$col[3],B=soton$col[5],C=soton$col[11])
  for(track in c("C","B","A")){
    ymin = ymin + onex/8
    text(0,ymin+onex/16,track,pos=2,xpd=TRUE,cex=0.8,col=trackcol[[track]])
    for(hapacc in hapdata[hapdata$Locus==locus & hapdata$Track==track,]$HapAcc){
      #trackdepdata = hapdepdata[hapdepdata$Locus==hapacc,]
      trackdepdata = hapdepdata[hapdepdata$HapAcc==hapacc,]
      trackhapdata = hapdata[hapdata$HapAcc==hapacc,]
      regionXPlot2(trackdepdata,trackhapdata,onex,ymin,trackcol[[track]])
    }
  }
  lines(plotdata$Pos/scaling,plotdata$X,col="blue",type="l",lwd=3)
  snpdata = hapsnpdata[hapsnpdata$Locus == locus,]
  if(dim(snpdata)[1] > 0){
    for(i in 1:dim(snpdata)[1]){
      # Locus	Block	Pos	Ref	A	B	nA	nB	pA
      snp = snpdata[i,]
      x = c(snp$Pos/scaling,snp$Pos/scaling)
      if(as.character(snp$Ref) == as.character(snp$A)){
        lines(x,c(0,snp$nB),type="l",col=trackcol$B)
      }
      if(as.character(snp$Ref) == as.character(snp$B)){
        lines(x,c(0,snp$nA),type="l",col=trackcol$A)
      }
      if(as.character(snp$Ref) != as.character(snp$A) & as.character(snp$Ref) != as.character(snp$B)){
        lines(x,c(0,(snp$nA+snp$nB)),type="l",col=trackcol$C)
      }
    }
  }
  abline(h=onex,col="red")
}



################# ::: GENOME FEATURE PLOTTING FUNCTIONS ::: ######################
#i# Features should be loaded into data frame, which is stored in a ftdb list object
#i# The list object should have a text key, e.g. "Reference" or "Assembly"
#i# The important fields for the input file are:
  #># locus = chromosome containing feature
  #># feature = feature type
  #># start, end = start and end positions 
#i# Other possible (unused) fields include: position, product, gene_synonym, note, db_xref, locus_tag, details
#i# At some point, a feature tag (locus_tag? name?) might also be used if present.

#i# The featureDefaults() function should be called prior to commandline arguments updating settings
featureDefaults = function(settings){
  settings$features = TRUE    # Whether to expect Features table plot.
  # Set list of features to include in plots
  settings$plotft = c("gene","mRNA","CDS","rRNA","tRNA","ncRNA","mobile","LTR","origin","centromere","telomere")
  return(settings)
}  

#i# The featureColours() function should be called AFTER commandline arguments updating settings
featureColours = function(settings){
  # Set up colours list if required
  if(is.null(settings$col)){
    settings$col = list()   # This will have contigs and chromosomes added elsewhere. 
  }
  # Add feature colours to list
  for(fi in 1:length(settings$plotft)){
    settings$col[[settings$plotft[fi]]] = soton$col[fi+1]
  }
  # Specific settings for LTRs and mobile elements. (Consider doing this for other things too!)
  settings$col[["LTR"]] = soton$col[15]
  settings$col[["mobile"]] = soton$col[16]
  return(settings)
}

#i# Want to have a generic function that loads features.
#i# Global ftdb list must already exist. (Cannot make global object within function.)
#># refbase = base of Feature.tdt file from which features to be loaded
#># ftype ["Reference"] = key for ftdb table
featureSetup = function(refbase,ftype="Reference"){
  # Load features
  if(is.null(settings$features)){ settings$features = TRUE }
  if("ftdb" %in% ls() == FALSE){ assign("ftdb", list(),  envir = .GlobalEnv) }
  ftfile = paste(refbase,".Feature.tdt",sep="")
  if(settings$features && file.exists(ftfile)){
    writeLines(paste("Adding",ftype,"to ftdb"))
    ftdata = read.table( ftfile, header = T, stringsAsFactors = F, sep = '\t', quote = '', comment.char="")
    ftdata = ftdata[ftdata$feature != "source",]
    ftdata = ftdata[ftdata$feature != "misc_RNA",]
    ftdata = ftdata[ftdata$feature != "misc_feature",]
    if(length(ftdata[ftdata$feature == "mobile_element",]$feature)>0){
      ftdata[ftdata$feature == "mobile_element",]$feature = "mobile"    
    }
    if(length(ftdata[ftdata$feature == "rep_origin",]$feature)>0){
      ftdata[ftdata$feature == "rep_origin",]$feature = "origin"
    }
    ftdata$feature = as.factor(ftdata$feature)
    colnames(ftdata)
    summary(ftdata)
    (levels(ftdata$feature))
  }else{ ftdata = data.frame() }
  ftdb[[ftype]] = ftdata 
  assign("ftdb", ftdb,  envir = .GlobalEnv)
  return(ftdata)
}
#i# Simple function that returns whether a given ftype has any features data loaded
#i# Use: if(isFeatureData(ftype)){ ...
isFeatureData = function(ftype){
  if(is.null(settings$features)){ settings$features = TRUE }
  if(settings$features == FALSE){ return(FALSE) }
  #!# #?# Why does the following always give FALSE even when it should be TRUE?! #!#
  #if("ftdb" %in% ls() == FALSE){ return(FALSE) }
  if(ftype %in% names(ftdb)){ 
    if(dim(ftdb[[ftype]])[1] > 0){
      return(TRUE)
    }
  }
  return(FALSE)
}

#i# Feature plot, plots features within xmin-xmax window for given chromosome
#># ftype = key of ftdb list
#># locus = locus name (matching locus in ftdb[[ftype]]
#># xmin, xmax = sets position limits of the plots, already scaled (e.g. kb)
#># acconly [TRUE] = whether to extract just the accession number from gene_SPEC__AccNum name
featurePlot = function(ftype,locus,xmin,xmax,acconly=TRUE){
  # Plot of the Features on a Chromosome:
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  plotft = settings$plotft
  ymax = length(plotft)
  features = ftdb[[ftype]]
  if(acconly){ locus = strsplit(locus,"[_.]")[[1]][4] }
  # Setup Plot
  plot(c(xmin,xmax),c(0,ymax),col="red",type="n",main=paste(locus,"features"),ylab="",xlab=paste("Chromosome Position",units),xaxs="i",yaxt="n")
  plotGridlines(xmin,xmax,0,ymax,10000/scaling,50000/scaling,1,ymax)
  # Plot Features
  if(settings$features){
    for(fi in 1:ymax){
      y = ymax-fi
      ftype = plotft[fi]
      text(xmin,y+0.5,ftype,pos=2,xpd=TRUE,cex=0.8)
      ftdata = features[features$feature == ftype & features$locus == locus,]
      if(dim(ftdata)[1] > 0){
        rect(ftdata$start/scaling,y,ftdata$end/scaling,y+1,col=settings$col[[ftype]],border=NA)  
      }
    }
  }
}




