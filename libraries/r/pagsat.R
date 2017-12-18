########################################################
### Pairwise Assembled Genome Sequence Analysis Tool ###
### VERSION: 2.0.1                             ~~~~~ ###
### LAST EDIT: 03/10/16                        ~~~~~ ###
### AUTHORS: Richard Edwards 2016              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

################# ::: HISTORY ::: ######################
# v0.4.0 : Initial version based on PacBio.R v0.3.0.
# v0.5.0 : Updated for PAGSAT v1.2.0.
# v0.6.0 : Altered identity to be relative to coverage, not length
# v0.7.0 : Added ChromAlignLoc output
# v0.7.1 : Added haploid/doploid modes for ChromAlignLoc output.
# v0.7.2 : Switched off some outputs if tables missing.
# v0.7.3 : Corrected crash in absence of Bwd hits.
# v0.7.4 : Added comment.char="" to read.table function calls to avoid issue with # in input. (Odd R default!!)
# v1.0.0 : Updating documentation in parallel with PAGSAT v1.12.0. Renamed to go with PAGSAT_V1
# NOTE: PAGSAT_V1 will now call frozen pagsat_V1.R file.
# v2.0.0 : Major overhaul in line with PAGSAT v2.0.0.
# v2.0.1 : Fixed a bug where there are no hits for a chromosome!

############### ::: GENERAL SETUP ::: ##################
# Usage = Rscript rje.r pagsat basefile [options] 
# If testing, make sure rje.r has been run to set up general colours (rje_col) and functions (rje_misc)

setTesting = function(){
  arglist = list("rtype"="pacbio","basefile"="/Users/redwards/Data/projects/YeastPacBio-Jun15/analysis/2015-08-14.PAGSAT.MBG8150/MBG8150.SP16496.hcq.qv20.sgd.srt.PAGSAT/MBG8150.SP16496.hcq.qv20.sgd.srt.L1000ID99","refbase"="/Users/redwards/Documents/Projects/YeastPacBio-Jun15/data/2015-07-01.S288C.8150.Ref/sgd.srt")
  arglist$basefile = "/Users/redwards/Data/projects/YeastPacBio-Jun15/dev/2016-09-26.PAGSATV1.12/MBG482001.sgd.ref.PAGSAT/MBG482001.sgd.ref.L250ID95"
  arglist$refbase = "/Users/redwards/Data/projects/YeastPacBio-Jun15/data/2016-01-19.Ref.S288C/sgd.ref"
  arglist$basefile = "/Users/redwards/Documents/Projects/YeastPacBio-Jun15/dev/2016-09-26.PAGSATV1.12/MBG482001.sgd.ref.PAGSAT/MBG482001.sgd.ref.L250ID95"
  arglist$refbase = "/Users/redwards/Documents/Projects/Microbiogen-Nov14/paper1/2016-01-19.Ref.S288C/sgd.ref"
  arglist$assbase = "/Users/redwards/Documents/Projects/YeastPacBio-Jun15/analysis/2016-09-22.MBG482001/MBG482001"
}

#i# The main setup of input data has now been moved to pagsat-setup.R
#i# arglist should already have been setup. (Either from rje.r or the setTesting() method above.)
rjesource("pagsat-setup.R")

#># seqdb has Sequences.tdt summary data
#># |-- rseqdb and aseqdb are Reference and Assembly subsets, respectively.
#># summdb has Summary Coverage data

############################### ::: DEFINE METHODS ::: ###########################
settings$testing = FALSE   # If true, will plot figure after each method

### ~ Call rje_genomics function definitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# This has some general chromosome plot functions
rjesource("rje_genomics.r")   

rjesource("pagsat-methods.R")

############# ::: MAIN RUN CODE ::: ##################


### ~ Generate multipanel plots for each reference chromosome ~~~~~~~~~~~~~ ###
for(chrname in levels(covplotdb[["Reference"]]$Qry)){
  pngfile=paste(settings$pngbase,"covplot",chrname,"png",sep=".")
  #CairoPNG(filename=pngfile, width=pwidth, height=pheight, pointsize=12)
  png(filename=pngfile, width=settings$pngwidth, height=settings$pngheight, units = "px", pointsize=settings$pointsize)
  panels = c(rep(1,8),rep(2,5),rep(3,5),rep(4,5))
  layout(matrix(panels,byrow=TRUE,nrow=length(panels)))
  par(mar=c(2,6,2,1))
  chromPlot(covplotdb[["Reference"]],chrname,minlablen=settings$minloclen)
  featurePlot(covplotdb[["Reference"]],ftdb,chrname)
  # Generate hitlist of contigs that have unique coverage 
  chrid = strsplit(chrname,"_")[[1]][1]
  mapcontigs = mapdb[mapdb$Chrom==chrid,]
  (hitlist = as.character(mapcontigs$Qry[order(mapcontigs$Unique,decreasing=TRUE)]))
  hitx = min(length(hitlist),settings$topcontigs)
  if(hitx > 0){
    localChromPlot(locdb,chrname,hitlist[1:hitx],ctype=paste("top",settings$topcontigs))
  }else{
    localChromPlot(locdb,chrname,hitlist,ctype=paste("top",settings$topcontigs))
  }
  par(mar=c(5,6,2,1))
  chromDifPlot(covplotdb[["Reference"]],chrname)
  dev.off()
}

### ~ Generate multipanel assembly plots for each reference chromosome ~~~~~~~~~~~~~ ###
for(chrname in levels(covplotdb[["Reference"]]$Qry)){
  # Setup data for plot. Need hitlist to determine size of plot and bottom panel
  chrid = strsplit(chrname,"_")[[1]][1]
  mapcontigs = mapdb[mapdb$Chrom==chrid,]
  # Generate hitlist of contigs that have unique coverage 
  (hitlist = as.character(mapcontigs$Qry[order(mapcontigs$Unique,decreasing=TRUE)]))
  # Setup PNG file  
  pheight = settings$pngheight
  hitrep = max(8,as.integer(length(hitlist)/2))
  #i# 6+5+3 = 14 versus 10 is the "standard" ratio
  topheight = as.integer(pheight * 14/24)
  botheight = as.integer(pheight * hitrep/24)
  pheight = topheight + botheight
  pngfile=paste(settings$pngbase,"assembly",chrname,"png",sep=".")
  #CairoPNG(filename=pngfile, width=pwidth, height=pheight, pointsize=12)
  png(filename=pngfile, width=settings$pngwidth, height=pheight, units = "px", pointsize=settings$pointsize)
  panels = c(rep(1,6),rep(2,5),rep(3,3),rep(4,hitrep))  # Entire bottom panel to have room for contigs
  layout(matrix(panels,byrow=TRUE,nrow=length(panels)))
  par(mar=c(2,6,2,1))
  #?#chromPlot(cnvdb,chrname,minlablen=settings$minloclen,recdata=FALSE)
  if(dim(covplotdb[["RefSnap"]])[1] > 0){
    chromPlot(covplotdb[["RefSnap"]],chrname,minlablen=settings$minloclen,recdata=FALSE)
  }else{
    chromPlot(covplotdb[["Reference"]],chrname,minlablen=settings$minloclen,recdata=FALSE)
  }
  featurePlot(covplotdb[["Reference"]],ftdb,chrname)
  chromosomeMap(rseqdb,uniqdb,ftdb,chromlist=c(chrname),pureplot = TRUE)
  par(mar=c(5,6,2,1))
  localChromPlot(uniqdb,chrname,hitlist,ctype=paste("top",settings$topcontigs))
  dev.off()
}



### ~ Generate multipanel assembly plots for each Assembly contig ~~~~~~~~~~~~~ ###
(medianx = median(xcovdb[["Assembly"]]$MedianX))
for(chrname in levels(covplotdb[["Assembly"]]$Qry)){
  pngfile=paste(settings$pngbase,"assembly",chrname,"png",sep=".")
  #CairoPNG(filename=pngfile, width=pwidth, height=pheight, pointsize=12)
  png(filename=pngfile, width=settings$pngwidth, height=settings$pngheight, units = "px", pointsize=settings$pointsize)
  #panels = c(rep(1,6),rep(2,5),rep(3,3),rep(4,10))  # Entire bottom panel to have room for contigs
  #1# Normal hit plot
  #2# Contig Depth Plot
  #3# Hitlist depth plots
  panels = c(rep(1,5),rep(2,5),rep(3,10))  # Entire bottom panel to have room for contigs
  layout(matrix(panels,byrow=TRUE,nrow=length(panels)))
  par(mar=c(2,6,2,1))

  #1# Normal hit plot
  chromPlot(covplotdb[["Assembly"]],chrname,minlablen=settings$minloclen,recdata=FALSE)

  #2# Contig Depth Plot
  if(settings$depth){
    chromDepthPlot(depthdb[["Assembly"]],chrname,medianx)
  }else{
    chromDifPlot(covplotdb[["Assembly"]],chrname)    
  }
  
  #!# At some point, add a plot of gene syntenty as features
  #!#featurePlot(covplotdb[["Reference"]],ftdb,chrname)
  # Generate hitlist of contigs that have unique coverage to same chromosome
  (chrhits = as.character(mapdb[mapdb$Qry==chrname,]$Chrom))
  maphits = c()
  for(chrom in chrhits){
    maphits = union(maphits,as.character(mapdb[mapdb$Chrom==chrom,]$Qry))
  }
  (maphits = maphits[maphits != chrname])

  locdata = asslocdb  #-># locdb[["RecAss"]]
  covsum = c()
  hitlist = c()
  if(length(maphits)>0){
    for(contig in maphits){
      sdata = locdata[locdata$Qry==chrname & locdata$Hit==contig & locdata$Dirn=="Fwd",]
      covsum = c(covsum,sum(sdata$Length))
    }
    (covsum)
    (covlist = maphits[order(covsum,decreasing=TRUE)])
    (covsum = covsum[order(covsum,decreasing=TRUE)])
    (hitlist = covlist[covsum>0])
    if(length(hitlist)>settings$topcontigs){
      hitlist = hitlist[1:settings$topcontigs]
    }
    
  }
  
  par(mar=c(5,6,2,1))
  if(length(hitlist) > 0){
    if(settings$depth){
      locAlnDepthPlot(chrname,hitlist,locdata,ctype="Assembly")
    }else{
      localChromPlot(locdata,chrname,hitlist,ctype="Assembly")
    }
  }
  
  dev.off()
}



### ~ Generate multipanel plots of top gene/protein hits for each reference chromosome ~~~~~~~~~~~~~ ###
if(settings$genesummary){
  for(chrname in levels(covplotdb[["Reference"]]$Qry)){
    pngfile=paste(settings$pngbase,"genehits",chrname,"png",sep=".")
    #CairoPNG(filename=pngfile, width=pwidth, height=pheight, pointsize=12)
    png(filename=pngfile, width=settings$pngwidth, height=settings$pngheight, units = "px", pointsize=settings$pointsize)
    panels = c(rep(1,7),rep(2,4),rep(3,6))
    layout(matrix(panels,byrow=TRUE,nrow=length(panels)))
    par(mar=c(2,6,2,1))
    chromGenesPlot(covplotdb[["Reference"]],genedb,protdb,chrname,ctype="Chromosome")
    featurePlot(covplotdb[["Reference"]],ftdb,chrname)
    chrid = strsplit(chrname,"_")[[1]][1]
    hitlist = levels(as.factor(as.vector(genedb[genedb$Chrom==chrid & genedb$Synteny != "Orphan",]$Hit)))
    par(mar=c(5,6,2,1))
    localChromPlot(locdb,chrname,hitlist,ctype="GeneHit")
    dev.off()
  }
}
  
### ~ Generate multipanel plots for each assembled contig ~~~~~~~~~~~~~ ###
for(contig in levels(covplotdb[["Assembly"]]$Qry)){
  pngfile=paste(settings$pngbase,"covplot",contig,"png",sep=".")
  #CairoPNG(filename=pngfile, width=pwidth, height=pheight, pointsize=12)
  png(filename=pngfile, width=settings$pngwidth, height=settings$pngheight, units = "px", pointsize=settings$pointsize)
  par(mar=c(2,6,2,1))
  #panels = c(1,1,1,2,2,2,3,3,1,1,1,2,2,2,4,4,1,1,1,2,2,2,5,5,1,1,1,2,2,2,6,6)
  #layout(matrix(panels,byrow=FALSE,nrow=8))
  panels = c(1,2,3,1,2,4,1,2,5,1,2,6)
  layout(matrix(panels,byrow=FALSE,nrow=3))
  chromPlot(covplotdb[["Assembly"]],contig,ctype="Contig")
  localContigPlot(covplotdb[["Assembly"]],locdb,contig)
  # Generate dotplots versus best matching chromosomes.
  #!# Might want to re-output and use the dircov table from PAGSAT to sort by Unique coverage, not total?
  #i# By sorting rest by total coverage, get good sense of repeats but not necessarily translocations.
  covsum = c()
  for(chrname in levels(locdb$Qry)){
    sdata = locdb[locdb$Qry==chrname & locdb$Hit==contig,]
    covsum = c(covsum,sum(sdata$Length))
  }
  (covsum)
  (covlist = levels(locdb$Qry)[order(covsum,decreasing=TRUE)])
  if(length(covlist)>4){ covlist = covlist[1:4] }

  # Make sure mapped csome is in here
  mapcontigs = mapdb[mapdb$Qry==contig,]
  (chridlist = as.character(mapcontigs$Chrom[order(mapcontigs$Unique,decreasing=TRUE)]))
  hitlist = c()
  for(chrid in chridlist){ hitlist = c(hitlist,seqData(chrid,"Name","Chrom","Reference")) }
  (hitlist)
      
  if(length(hitlist)>4){ hitlist = hitlist[1:4] }
  if(length(hitlist)<1){ hitlist = covlist[1:4] }
  i = 0
  while(length(hitlist)<length(covlist)){
    i = i + 1
    addme = TRUE
    for(j in 1:length(hitlist)){
      if(hitlist[j] == covlist[i]){ addme = FALSE }
    }
    if(addme){ hitlist = c(hitlist,covlist[i]) } 
  }
  
  par(mar=c(5,6,1,1))
  for(chrname in hitlist){  # Have chrom as Query otherwise need to reverse Qry/Hit(Sbj) data - also allows alignments of contig
    if(dim(locdb[locdb$Qry==chrname & locdb$Hit==contig,])[1] > 0){
      dotplot(locdb,chrname,contig,qlen=max(covplotdb[["Reference"]][covplotdb[["Reference"]]$Qry==chrname,]$Pos),hlen=max(covplotdb[["Assembly"]][covplotdb[["Assembly"]]$Qry==contig,]$Pos),pcol=pcol,minfrag=settings$minloclen)      
    }else{
      plot(c(1,2),c(1,2),type="n",xlab="-",ylab="-",axes=TRUE,ann=TRUE,mar=c(0,1,4,1))      
    }    
  }
  dev.off()
}

############# ::: SUMMARY DATA CODE ::: ##################

## ~ SUMMARY DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## ~ Two-panel plot of overal summary and chromosome coverage and contig coverage ~ ##
pngfile=paste(settings$pngbase,"summary","png",sep=".")
#CairoPNG(filename=pngfile, width=pwidth, height=pheight, pointsize=12)
png(filename=pngfile, width=settings$pngwidth, height=settings$pngheight, units = "px", pointsize=settings$pointsize)
panels = c(rep(1,3),rep(2,5))
layout(matrix(panels,byrow=TRUE,nrow=length(panels)))
par(mar=c(5,6,1,1))
summGraph(summdb)
chromosomeMap(rseqdb,locdb,ftdb,gtype="Reference")
dev.off()


## ~ Plot contig coverage ~ ##
apngheight = settings$pngheight
apointsize = settings$pointsize
cnum = dim(aseqdb)[1]
if(as.integer(cnum/25) > as.integer(cnum/30)){
  cnum = min(dim(aseqdb)[1],30)
}else{
  cnum = min(dim(aseqdb)[1],25)
}
if(cnum > 50){
  apngheight = 2 * apngheight 
  apointsize = 0.5 * apointsize * cnum / 50
}else{
  if(cnum > 25){
    apngheight = cnum * apngheight / 25  
    apointsize = apointsize * 0.8
  }
}

## Will need multiple plots if too many contigs
if(cnum < dim(aseqdb)[1]){
  nplot = as.integer(dim(aseqdb)[1]/cnum) + 1
  (chrnames = aseqdb[order(aseqdb$Length,decreasing=TRUE),]$Qry)
  for(i in 1:nplot){
    cx = 1 + (cnum * (i - 1))
    cy = min(cnum * i,dim(aseqdb)[1])
    pngfile=paste(settings$pngbase,"assembly",i,"png",sep=".")
    png(filename=pngfile, width=settings$pngwidth, height=apngheight, units = "px", pointsize=apointsize)
    par(mar=c(5,6,1,1))
    chromosomeMap(aseqdb,locdb,ftdb,gtype="Assembly",chromlist=chrnames[cx:cy])
    dev.off()
  }    
}else{  
  pngfile=paste(settings$pngbase,"assembly","png",sep=".")
  png(filename=pngfile, width=settings$pngwidth, height=apngheight, units = "px", pointsize=apointsize)
  par(mar=c(5,6,1,1))
  chromosomeMap(aseqdb,locdb,ftdb,gtype="Assembly")
  dev.off()
}
