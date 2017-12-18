########################################################
### RJE Genomics Plot Functions                ~~~~~ ###
### VERSION: 0.2.0                             ~~~~~ ###
### LAST EDIT: 28/03/17                        ~~~~~ ###
### AUTHORS: Richard Edwards 2016              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for genomic plotting functions to be used by other scripts, 
# e.g. PAGSAT.R and samtools.R

################# ::: HISTORY ::: ######################
# v0.0.0 : Initial version based on PAGSAT v0.7.4 for use with samtools.R v0.0.0.
# v0.1.0 : Added SAMPhaser methods for used with SAMPhaser v0.4.0.
# v0.2.0 : Added Genome feature setup and plotting methods


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
samphaserStackPlot = function(depdata,hapdepdata,hapdata,locus,onex){
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

samphaserStackPlot2 = function(depdata,hapdepdata,hapdata,hapsnpdata,locus,onex,zoom=c(0,0)){
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
      trackdepdata = hapdepdata[hapdepdata$Locus==hapacc,]
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




