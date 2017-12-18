########################################################
### Pairwise Assembled Genome Sequence Analysis Tool ###
### VERSION: 2.0.1                             ~~~~~ ###
### LAST EDIT: 02/12/16                        ~~~~~ ###
### AUTHORS: Richard Edwards 2016              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

################# ::: HISTORY ::: ######################
# v2.0.0 : Major overhaul as part of PAGSAT v2.0.0.
# v2.0.1 : Fixed a bug where there are no hits for a chromosome!

############################### ::: DEFINE METHODS ::: ###########################

### ~ Main Summary Table plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#># summdata = Summary.tdt data frame
summGraph = function(summdata){
  # Setup Plot
  scaling = settings$scaling   # Used for units of text display
  xmax = 100                   # Printing as a percentage of total
  # sumtypes will cross-reference to summdata$Summary; sumtext will be printed
  sumtypes = c("Reference","Assembly","RefBench","Genes","Proteins","ReferenceAlign","AssemblyAlign")
  sumtext = c("Reference","Assembly","RefBench","Ref.Genes","Ref.Proteins","Ref.Align","AssemblyAlign")
  # Determine those to output. Will leave space for missing gene/protein data.
  ymax = 5   
  if(settings$chromalign){ ymax = 7 }
  # Generate plot area
  par(mar=c(5,6,2,1))
  plot(c(0,xmax),c(0,ymax),col="red",type="n",main=paste(settings$basename,"summary"),ylab="",xlab=paste("Percentage"),xaxs="i",yaxt="n")
  plotGridlines(0,xmax,0,ymax,5,10,0.5,1)
  # Plot Summary
  for(fi in 1:ymax){
    y = ymax-fi
    text(0,y+0.5,sumtext[fi],pos=2,xpd=TRUE,cex=0.8)
    ydata = summdata[summdata$Summary==sumtypes[fi],]
    if(dim(ydata)[1] > 0){
      rect(0,y+0.5,100.0*ydata$Coverage/ydata$Length,y+1,col=soton$col[12])
      text(0,y+0.75,paste(format(100.0*ydata$Coverage/ydata$Length,digits=5),"% coverage / ",(ydata$Length-ydata$Coverage)/scaling,settings$units," missing",sep=""),pos=4,xpd=TRUE,cex=0.8)
      rect(0,y+0.5,100.0*ydata$Identity/ydata$Coverage,y,col=soton$col[20])
      if((ydata$Coverage-ydata$Identity)>9999){
        text(0,y+0.25,paste(format(100.0*ydata$Identity/ydata$Coverage,digits=5),"% identity / ",format((ydata$Coverage-ydata$Identity)/1000,digits=3),"kb variant (",ydata$Perfect," of ",ydata$N," identical)",sep=""),pos=4,xpd=TRUE,cex=0.8)
      }else{          
        text(0,y+0.25,paste(format(100.0*ydata$Identity/ydata$Coverage,digits=5),"% identity / ",(ydata$Coverage-ydata$Identity),"bp variant (",ydata$Perfect," of ",ydata$N," identical)",sep=""),pos=4,xpd=TRUE,cex=0.8)
      }
    }else{
      text(0,y+0.5,"No data generated.",pos=4,xpd=TRUE,cex=0.8)
    }
  }
}
if(settings$testing){ summGraph(summdb) }
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###



### ~ Main Summary Chromosome Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i## Mapping of locdb onto all chromosomes, by size.
#i# Used to be called coverageByChromosome but now chromosomeMap
#># covdata = Coverage data to use (subset of Sequences.tdt data frame)
#># locdata = Local BLAST hit table (Reference vs Assembly)
#># features = Reference Feature Table
#># gtype="Reference" or "Assembly"
#># pureplot:bool [FALSE] = Whether to exclude the identity/coverage information and make the axis purely the length of the chromosome.
#i# pureplot=TRUE is designed for a single chromosome.
chromosomeMap = function(covdata,locdata,features,gtype="Reference",chromlist=c(),pureplot=FALSE){
  # Setup plot region and scaling parameters 
  scaling = settings$scaling
  # Setup Empty Plot Region and Labels. 
  if(pureplot){
    xmax = max(covdata[covdata$Qry==chromlist,]$Length)/scaling
  }else{ 
    par(mar=c(5,6,2,1))
    xmax=max(covdata$Length)/scaling 
  }
  chrnames = covdata[order(covdata$Length,decreasing=TRUE),]$Qry
  if(length(chromlist) > 0){ chrnames = chromlist }
  ymax=length(chrnames)
  if(pureplot){
    plot(c(0,xmax),c(0,ymax),col="red",type="n",main=paste(settings$basename,gtype,"coverage"),ylab="",xlab=paste("Length (",settings$units,")",sep=""),xaxs="i",yaxs="i",yaxt="n")
  }else{
    plot(c(0,xmax*1.2),c(0,ymax),col="red",type="n",main=paste(settings$basename,gtype,"coverage"),ylab="",xlab=paste("Length (",settings$units,")",sep=""),xaxs="i",yaxs="i",yaxt="n")
  }
  plotGridlines(0,xmax,0,ymax,50,200,0.25,1)
  # Setup list of chromosomes, abbrev csome names and loci for features ("Reference" only)
  csomes = c()
  clocus = c()
  for(chrname in chrnames){
    csomes = c(csomes,strsplit(chrname,"_")[[1]][1])
    clocus = c(clocus,strsplit(chrname,"[_.]")[[1]][4])
  }
  # Plot chromosomes and, features and local mappings
  for(fi in 1:ymax){
    y = ymax-fi
    xlen = covdata[covdata$Qry==chrnames[fi],"Length"]/scaling
    #chromcol = settings$col[[as.character(chrnames[fi])]]
    chromcol = chromCol(chrnames[fi])
    ydata = covdata[covdata$Qry==chrnames[fi],]
    text(0,y+0.5,csomes[fi],pos=2,xpd=TRUE,cex=0.8,col=chromcol)
    if(pureplot == FALSE){
      text(xlen,y+0.75,paste(format(100.0*ydata$Coverage/ydata$Length,digits=5),"% coverage / ",(ydata$Length-ydata$Coverage)/scaling,settings$units," missing",sep=""),pos=4,xpd=TRUE,cex=max(0.4,1 - 0.015*ymax))
      text(xlen,y+0.25,paste(format(100.0*ydata$Identity/ydata$Coverage,digits=5),"% identity / ",(ydata$Coverage-ydata$Identity),"bp variant",sep=""),pos=4,xpd=TRUE,cex=max(0.4,1 - 0.015*ymax))
    }
    # Plot chromosome length and (possibly) feature data
    if(gtype == "Reference"){
      ldata = locdata[locdata$Qry == as.character(chrnames[fi]),]
    }else{
      ldata = locdata[locdata$Hit == as.character(chrnames[fi]),]
    }
    if(dim(features[features$locus == clocus[fi],])[1] > 0){
      rect(0,y+0.5,xlen,y+0.75,col="white")
      for(ftype in settings$plotft[settings$plotft != "mRNA"]){
        ftdata = features[features$feature == ftype & features$locus == clocus[fi],]
        if(dim(ftdata)[1] > 0){
          rect(ftdata$start/scaling,y+0.5,ftdata$end/scaling,y+0.75,col=settings$col[[ftype]],border=NA)  
        }
      }      
    }else{
      rect(0,y+0.5,xlen,y+0.75,col=chromcol)
    }
    # Plot mapped local hits
    ldata = ldata[order(ldata$Length),]
    #i# Need to use $HitCol and $QryCol here because each rect might be a different chromosome/contig
    if(gtype == "Reference"){
      for(chrname in seqdb[seqdb$Source=="Assembly","Name"]){
        dirndata = ldata[ldata$Dirn == "Fwd" & ldata$Hit == chrname,]
        if(dim(dirndata)[1] > 0){
          rect(dirndata$QryStart/scaling,y+0.75,dirndata$QryEnd/scaling,y+1,col=chromCol(chrname)) 
        }
        dirndata = ldata[ldata$Dirn == "Bwd" & ldata$Hit == chrname,]
        if(dim(dirndata)[1] > 0){
          rect(dirndata$QryStart/scaling,y+0.25,dirndata$QryEnd/scaling,y+0.5,col=chromCol(chrname))  
        }
      }
    }else{
      for(chrname in seqdb[seqdb$Source=="Reference","Name"]){
        dirndata = ldata[ldata$Dirn == "Fwd" & ldata$Qry == chrname,]
        if(dim(dirndata)[1] > 0){
          rect(dirndata$SbjStart/scaling,y+0.75,dirndata$SbjEnd/scaling,y+1,col=chromCol(chrname))
        }
        dirndata = ldata[ldata$Dirn == "Bwd" & ldata$Qry == chrname,]
        if(dim(dirndata)[1] > 0){
          rect(dirndata$SbjStart/scaling,y+0.25,dirndata$SbjEnd/scaling,y+0.5,col=chromCol(chrname))  
        }
      }
    }
  }
}
if(settings$testing){ chromosomeMap(rseqdb,locdb,ftdb) }
if(settings$testing){ chromosomeMap(rseqdb,uniqdb,ftdb) }
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###













## ~ Chromosome gene/protein TopHits plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
chromGenesPlot = function(covplotdata,genedata,protdata,chrname,ctype="Chromosome"){
  # Plot of the Gene/Protein Hits betweeen Chromosomes and Contigs
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  plotdata = covplotdata[covplotdata$Qry==chrname,]
  xmax = max(plotdata$Pos)/scaling
  ymax = 100.0
  ymin = 0
  chrid = strsplit(chrname,"_")[[1]][1]
  geneplot = genedata[genedata$Qry==chrid,]  # If plotting a contig, will want to change this to $Hit
  geneplot$x = (geneplot$Start+geneplot$End)/2
  protplot = protdata[protdata$Qry==chrid,]  # If plotting a contig, will want to change this to $Hit
  protplot$x = (protplot$Start+protplot$End)/2
  # Setup Plot
  plot(c(0,xmax),c(ymin,ymax),col="red",type="n",main=paste(settings$basename,"-",chrname,"-",xmax,settings$units),ylab="Top gene/protein hits (%identity)",xlab=paste(ctype,"Position",units),xaxs="i")
  plotGridlines(0,xmax,ymin,ymax,10000/scaling,50000/scaling,5,20)
  # Plot Data
  for(xcontig in levels(geneplot$Hit)){
    cdata = geneplot[geneplot$Hit==xcontig & geneplot$Synteny == "Orphan",]
    if(xcontig == ""){ 
      xcontig = "NULL" 
    }
    points(cdata$x/scaling,cdata$Qry_AlnID,col="black",pch=21,bg=settings$col[[xcontig]])
  }
  for(xcontig in levels(protplot$Hit)){
    cdata = protplot[protplot$Hit==xcontig & protplot$Synteny == "Orphan",]
    if(xcontig == ""){ 
      xcontig = "NULL" 
    }
    points(cdata$x/scaling,cdata$Qry_AlnID,col="black",pch=23,bg=settings$col[[xcontig]])
  }
  for(xcontig in levels(geneplot$Hit)){
    cdata = geneplot[geneplot$Hit==xcontig & geneplot$Synteny != "Orphan",]
    if(xcontig == ""){ 
      xcontig = "NULL" 
    }
    points(cdata$x/scaling,cdata$Qry_AlnID,col=settings$col[[xcontig]],pch=1)
  }
  for(xcontig in levels(protplot$Hit)){
    cdata = protplot[protplot$Hit==xcontig & protplot$Synteny != "Orphan",]
    if(xcontig == ""){ 
      xcontig = "NULL" 
    }
    points(cdata$x/scaling,cdata$Qry_AlnID,col=settings$col[[xcontig]],pch=0)
  }
  points(xmax/20,12.5,col="black",pch=0,cex=2)
  text(xmax/18,12.5,"Protein",pos=4)
  points(xmax/20,7.5,col="black",pch=1,cex=2)
  text(xmax/18,7.5,"Gene",pos=4)
}

## ~ Chromosome coverage plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
#!# Need to change ctype to be "Reference" or "Assembly"
#># recdata:bool [TRUE] = Whether to plot Reciprocal data as negative axis. Otherwise, just plot positive data.
chromPlot = function(covplotdata,chrname,ctype="Chromosome",minlablen=1000,recdata=TRUE){
  # Plot of the Hits betweeen Chromosomes and Contigs(+) or Chromosomes(-):
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  plotdata = covplotdata[covplotdata$Qry==chrname,]
  plotdata = plotdata[(order(plotdata$Pos)),]
  xmax = max(plotdata$Pos)/scaling
  if(recdata){ ymin = -max(plotdata$RecAln) }
  else{ ymin = 0 }
  ymax = max(plotdata$HitAln,5)
  # Setup Plot
  plot(c(0,xmax),c(ymin,ymax),col="red",type="n",main=paste(settings$basename,"-",chrname,"-",xmax,settings$units),ylab="BLAST hits",xlab=paste(ctype,"Position",units),xaxs="i")
  plotGridlines(0,xmax,ymin,ymax,10000/scaling,50000/scaling,5,20)
  # Plot Data
  lines(plotdata$Pos/scaling,plotdata$HitAln,col="red",type="l",lwd=2)
  lines(plotdata$Pos/scaling,plotdata$HitSeq,col="blue",type="l",lwd=2)
  if(recdata){
    lines(plotdata$Pos/scaling,-plotdata$RecAln,col="red",type="l",lwd=2)
    lines(plotdata$Pos/scaling,-plotdata$RecSeq,col="blue",type="l",lwd=2)
  }
  for(xi in 2:length(plotdata$Pos)){
    xcontig = as.character(plotdata$HitChrom[xi])
    if(plotdata$Class[xi-1] == "U" & plotdata$Class[xi] == "U"){
      rect(plotdata$Pos[xi-1]/scaling,-1,plotdata$Pos[xi]/scaling,1,col=settings$col[[xcontig]],border="blue",lwd=1)
      if((plotdata$Pos[xi] - plotdata$Pos[xi-1]) >= minlablen){
        text((plotdata$Pos[xi-1]+plotdata$Pos[xi])/(2*scaling),max(2,ymax/20),xcontig,pos=4,srt=90,offset=0,col=settings$col[[xcontig]],xpd=TRUE)
      }
    }
    if(plotdata$Class[xi-1] == "C" & plotdata$Class[xi] == "C"){
      rect(plotdata$Pos[xi-1]/scaling,-1,plotdata$Pos[xi]/scaling,1,col=settings$col[[xcontig]],border="blue",lwd=1)
      if((plotdata$Pos[xi] - plotdata$Pos[xi-1]) >= minlablen){
        text((plotdata$Pos[xi-1]+plotdata$Pos[xi])/(2*scaling),max(2,ymax/20),xcontig,pos=4,srt=90,offset=0,col=settings$col[[xcontig]],xpd=TRUE)
      }
    }
    if(plotdata$Class[xi-1] == "N" & plotdata$Class[xi] == "N"){
      points((plotdata$Pos[xi-1]+plotdata$Pos[xi])/(2*scaling),max(1,ymax/10),pch=2,col=soton$col[5],lwd=2)
      text((plotdata$Pos[xi-1]+plotdata$Pos[xi])/(2*scaling),max(2.5,ymax/2),paste((plotdata$Pos[xi]-plotdata$Pos[xi-1])/scaling,settings$units),col=soton$col[5],pos=4,srt=90,xpd=TRUE,offset=0)
    }
  }
  #Legend
  legy = 0.95 * ymin
  legx = xmax/20
  text(0,0.95 * ymax,"Assembly vs Reference",pos=4)
  if(recdata){
    if(ctype == "Chromosome"){
      text(0,legy,"Reference vs Reference",pos=4)
    }else{
      text(0,legy,"Assembly vs Assembly",pos=4)
    }
  }else{ legy = 0.95 * ymax }
  lines(c(8:9)*legx,c(legy,legy),col="red",type="l",lwd=2)
  text(9*legx,legy,"Local",pos=4)
  lines(c(10:11)*legx,c(legy,legy),col="blue",type="l",lwd=2)
  text(11*legx,legy,"Chromosome/Contig",pos=4)
}

## ~ Chromosome coverage difference plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
chromDifPlot = function(covplotdata,chrname){
  # Plot of the Hits betweeen Chromosomes and Contigs(+) or Chromosomes(-):
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  plotdata = covplotdata[covplotdata$Qry==chrname,]
  plotdata = plotdata[(order(plotdata$Pos)),]
  xmax = max(plotdata$Pos)/scaling
  ymin = min(plotdata$HitAln-plotdata$RecAln)
  ymax = max(plotdata$HitAln-plotdata$RecAln)
  # Setup Plot
  plot(c(0,xmax),c(ymin,ymax),col="red",type="n",ylab="BLAST Hits",main="Assembly - Reference",xlab=paste("Chromosome Position",units),xaxs="i")
  plotGridlines(0,xmax,ymin,ymax,10000/scaling,50000/scaling,2,10)
  # Plot Data
  lines(plotdata$Pos/scaling,plotdata$HitAln-plotdata$RecAln,col="red",type="l")
  lines(plotdata$Pos/scaling,plotdata$HitSeq-plotdata$RecSeq,col="blue",type="l")
}

## ~ Features plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
featurePlot = function(covplotdata,features,chrname){
  # Plot of the Features on a Chromosome:
  locus = strsplit(chrname,"[_.]")[[1]][4]
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  plotdata = covplotdata[covplotdata$Qry==chrname,]
  plotdata = plotdata[(order(plotdata$Pos)),]
  xmax = max(plotdata$Pos)/scaling
  plotft = c("gene","mRNA","CDS","rRNA","tRNA","ncRNA","mobile","LTR","origin","centromere","telomere")
  ymax = length(plotft)
  # Setup Plot
  plot(c(0,xmax),c(0,ymax),col="red",type="n",main=paste(chrname,"features"),ylab="",xlab=paste("Chromosome Position",units),xaxs="i",yaxt="n")
  plotGridlines(0,xmax,0,ymax,10000/scaling,50000/scaling,1,ymax)
  # Plot Features
  if(settings$features){
    for(fi in 1:ymax){
      y = ymax-fi
      ftype = plotft[fi]
      text(0,y+0.5,ftype,pos=2,xpd=TRUE,cex=0.8)
      ftdata = features[features$feature == ftype & features$locus == locus,]
      if(dim(ftdata)[1] > 0){
        rect(ftdata$start/scaling,y,ftdata$end/scaling,y+1,col=settings$col[[ftype]],border=NA)  
      }
    }
  }
}




## ~ Local stack plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
localChromPlot = function(locdata,chrname,hitlist,ctype="top"){
  # Plot of the Features on a Chromosome:
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  xmax = seqData(chrname,"Length")/scaling
  ymax = max(1,length(hitlist))
  # Setup Plot
  plot(c(0,xmax),c(0,ymax),col="red",type="n",main=paste(chrname,ctype,"contigs hit stack"),ylab="",xlab=paste("Chromosome Position",units),xaxs="i",yaxt="n")
  if(length(hitlist) < 1){
    text(xmax/2,ymax/2,"No Hits",xpd=TRUE,cex=0.8)
    return()
  }
  plotGridlines(0,xmax,0,ymax,10000/scaling,50000/scaling,1,3)
  # Plot Features
  for(fi in 1:ymax){
    y = ymax-fi
    hit = hitlist[fi]
    hitname = strsplit(hit,"_")[[1]][1]
    hitlen = seqData(hit,"Length")
    text(0,y+0.5,hitname,pos=2,xpd=TRUE,cex=0.8)
    #!# ftdata as a name is confusing! These are local alignment chunks -> alndata
    ftdata = locdata[locdata$Qry == chrname & locdata$Hit == hit,]
    ftdata = ftdata[order(ftdata$Length),]
    dirndata = ftdata[ftdata$Dirn == "Fwd",]
    if(dim(dirndata)[1] > 0){
      rect(dirndata$QryStart/scaling,y+0.5,dirndata$QryEnd/scaling,y+1,col=chromCol(hit))  
      enddata = dirndata[dirndata$SbjStart == 1,]
      points(enddata$QryStart/scaling,rep(y+1.1,dim(enddata)[1]),pch=23,bg=chromCol(hit))
      enddata = dirndata[dirndata$SbjEnd == hitlen,]
      points(enddata$QryEnd/scaling,rep(y+1.1,dim(enddata)[1]),pch=23,bg=chromCol(hit))
    }
    dirndata = ftdata[ftdata$Dirn == "Bwd",]
    if(dim(dirndata)[1] > 0){
      rect(dirndata$QryStart/scaling,y,dirndata$QryEnd/scaling,y+0.5,col=chromCol(hit))  
    }
  }
}


# This is the full refereence stack. Not sure that it needs covplot data at all!
localContigPlot = function(covplotdata,locdata,contig){
  # Plot of the Features on a Chromosome:
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  plotdata = covplotdata[covplotdata$Qry==contig,] 
  plotdata = plotdata[(order(plotdata$Pos)),]
  xmax = max(plotdata$Pos)/scaling
  ymax = length(levels(locdata$Qry))
  # Setup Plot
  plot(c(0,xmax),c(0,ymax),col="red",type="n",main=paste(contig,"hit stack"),ylab="",xlab=paste("Contig Position",units),xaxs="i",yaxt="n")
  plotGridlines(0,xmax,0,ymax,10000/scaling,50000/scaling,0.5,1)
  # Plot Features
  for(fi in 1:ymax){
    y = ymax-fi
    chrname = levels(locdata$Qry)[fi]
    chrid = strsplit(chrname,"_")[[1]][1]
    text(0,y+0.5,chrid,pos=2,xpd=TRUE,cex=0.8)
    ftdata = locdata[locdata$Hit == contig & locdata$Qry == chrname,]
    ftdata = ftdata[order(ftdata$Length),]
    dirndata = ftdata[ftdata$Dirn == "Fwd",]
    if(dim(dirndata)[1] > 0){
      rect(dirndata$SbjStart/scaling,y+0.5,dirndata$SbjEnd/scaling,y+1,col=settings$col[[chrname]])  
    }
    dirndata = ftdata[ftdata$Dirn == "Bwd",]
    if(dim(dirndata)[1] > 0){
      rect(dirndata$SbjStart/scaling,y,dirndata$SbjEnd/scaling,y+0.5,col=settings$col[[chrname]])  
    }
  }
}

## ~ DotPlot of local alignments ~~~~~~~~~~~~~~~~~~~~~~~~ ##
# 100% = black
# >= 99% = purple
# >= 90% = blue
# < 90% = red
(pcol = c(rep("red",89),rep("blue",9),"purple","black"))
dotplot = function(locdata,qry,hit,qlen,hlen,pcol=rep(soton$col[1],100),minfrag=0){
  ### ~ Setup Plot ~ ###
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  if(is.na(hlen) | is.na(qlen)){ return() }
  if(hlen <= 0 | qlen <= 0){ return() }
  #title(main=paste(hit," vs ",qry," (Minfrag=",minfrag,"bp)",sep=""))
  ### ~ Plot local hit data ~ ###
  ldata = locdata[locdata$Qry==qry & locdata$Hit==hit & locdata$Length >= minfrag,]
  ldata$Perc = ldata$Identity / ldata$Length
  qry = strsplit(qry,"_")[[1]][1]
  hit = strsplit(hit,"_")[[1]][1]
  plot(c(1,qlen)/scaling,c(1,hlen)/scaling,type="n",xlab=paste(qry,units),ylab=paste(hit,units),axes=TRUE,ann=TRUE,mar=c(0,1,4,1))
  plotGridlines(0,qlen/scaling,0,hlen/scaling,50000/scaling,200000/scaling,50000/scaling,200000/scaling,xcap=TRUE,ycap=TRUE)
  text(0,0.95*hlen/scaling,paste("minfrag =",minfrag),pos=4)
  for(i in 1:length(ldata$Perc)){
    (tran = as.integer(100*ldata$Perc[i]) - 1)
    #(ldata[i,])
    #icol = (paste(pcol,tran,sep=""))
    icol = pcol[tran+1]
    lines(c(ldata$QryStart[i],ldata$QryEnd[i])/scaling,c(ldata$SbjStart[i],ldata$SbjEnd[i])/scaling,type="l",lwd=2,col=icol)
    lines(c(ldata$QryStart[i],ldata$QryEnd[i])/scaling,c(ldata$SbjStart[i],ldata$SbjEnd[i])/scaling,type="p",lwd=2,col=icol)
  }
}



  


## ~ alignedChromosomes mapping of aligned locdb onto all chromosomes, by size ~ ##
alignedChromosomes = function(covdata,alndata,features){
  scaling = settings$scaling
  # Setup Plot
  xmax = max(covdata$Length)/scaling #- this will be for individual chromosomes
  ymax = dim(covdata)[1]
  par(mar=c(5,6,2,1))
  plot(c(0,xmax),c(0,ymax),col="red",type="n",main=paste(settings$basename,"aligned coverage"),ylab="",xlab=paste("Length (",settings$units,")",sep=""),xaxs="i",yaxt="n")
  plotGridlines(0,xmax,0,ymax,50,200,0.25,1)
  # Setup list of chromosomes, abbrev csome names and loci for features ("Reference" only)
  chrnames = covdata[order(covdata$Length,decreasing=TRUE),]$Qry
  csomes = c()
  clocus = c()
  for(chrname in chrnames){
    csomes = c(csomes,strsplit(chrname,"_")[[1]][1])
    clocus = c(clocus,strsplit(chrname,"[_.]")[[1]][4])
  }
  # Plot chromosomes and, features and local mappings
  for(fi in 1:ymax){
    y = ymax-fi
    xlen = covdata[covdata$Qry==chrnames[fi],"Length"]/scaling
    chromcol = settings$col[[as.character(chrnames[fi])]]
    ydata = covdata[covdata$Qry==chrnames[fi],]
    text(0,y+0.5,csomes[fi],pos=2,xpd=TRUE,cex=0.8,col=chromcol)
    #text(xlen,y+0.75,paste(format(100.0*ydata$Coverage/ydata$Length,digits=5),"% coverage / ",(ydata$Length-ydata$Coverage)/scaling,settings$units," missing",sep=""),pos=4,xpd=TRUE,cex=max(0.4,1 - 0.015*ymax))
    #text(xlen,y+0.25,paste(format(100.0*ydata$Identity/ydata$Length,digits=5),"% identity / ",(ydata$Coverage-ydata$Identity),"bp variant",sep=""),pos=4,xpd=TRUE,cex=1 - 0.015*ymax)
    #text(xlen,y+0.25,paste(format(100.0*ydata$Identity/ydata$Coverage,digits=5),"% identity / ",(ydata$Coverage-ydata$Identity),"bp variant",sep=""),pos=4,xpd=TRUE,cex=max(0.4,1 - 0.015*ymax))
    # Plot chromosome length and (possibly) feature data
    rect(0,y+0.5,xlen,y+0.75,col="white")
    for(ftype in settings$plotft[settings$plotft != "mRNA"]){
      ftdata = features[features$feature == ftype & features$locus == clocus[fi],]
      if(dim(ftdata)[1] > 0){
        rect(ftdata$start/scaling,y+0.5,ftdata$end/scaling,y+0.75,col=settings$col[[ftype]],border=NA)  
      }
    }   
  
    # Plot mapped local hits
    ldata = alndata[alndata$Query == as.character(chrnames[fi]),]
    ldata = ldata[order(ldata$Length),]
    # Cycle through Ctid and Dirn data
    dirndata = ldata[ldata$Dirn == "Fwd" & ldata$Ctid == "B",]
    if(dim(dirndata)[1] > 0){
      rect(dirndata$QryStart/scaling,y+0.875,dirndata$QryEnd/scaling,y+1,col=dirndata$HitCol)  
    }
    dirndata = ldata[ldata$Dirn == "Fwd" & ldata$Ctid != "B",]
    if(dim(dirndata)[1] > 0){
      if(settings$diploid){
        rect(dirndata$QryStart/scaling,y+0.75,dirndata$QryEnd/scaling,y+0.875,col=dirndata$HitCol)  
      }else{
        rect(dirndata$QryStart/scaling,y+0.75,dirndata$QryEnd/scaling,y+1,col=dirndata$HitCol)          
      }
    }
    dirndata = ldata[ldata$Dirn == "Bwd" & ldata$Ctid != "B",]
    if(dim(dirndata)[1] > 0){
      if(settings$diploid){
        rect(dirndata$QryStart/scaling,y+0.375,dirndata$QryEnd/scaling,y+0.5,col=dirndata$HitCol)  
      }else{
        rect(dirndata$QryStart/scaling,y+0.25,dirndata$QryEnd/scaling,y+0.5,col=dirndata$HitCol)  
      }
    }
    dirndata = ldata[ldata$Dirn == "Bwd" & ldata$Ctid == "B",]
    if(dim(dirndata)[1] > 0){
      rect(dirndata$QryStart/scaling,y+0.25,dirndata$QryEnd/scaling,y+0.375,col=dirndata$HitCol)  
    }
  }
}

#chromosomeMap(rseqdb,locdb,ftdb,gtype="Reference")
#chromosomeMap(aseqdb,locdb,ftdb,gtype="Assembly")


#!# Need a method to convert a local match into a coverage plot
#!# Map ends onto coverage data - if not in data then must be between two values
#!# -> Can just get all values in the Start-End range and then duplicate the end ones! (Easy!)

### ~ Methods to convert locdb coordinates from Qry to Hit or visa versa ~ ###
#?# Do we need this? Just load the RecAss local data?
#># ymin and ymax determine the vertical plot positions so that this code can be reused in different plots
#># ysplit sets the proportion of the yaxis that is the local plot, with the rest the xplot
locXPlot = function(qry,hit,locdata,depdata,maxdep,avdep,ymin,ymax,ysplit=0.2){
  scaling = settings$scaling
  # Set up local data to plot
  ldata = locdata[locdata$Qry == qry & locdata$Hit == hit & locdata$Dirn == "Fwd",]
  ldata = ldata[order(ldata$Length),]
  hitctg = seqData(hit,"Chrom")
  xdata = depdata[depdata$Locus == hit,]
  # Cycle through each alignment
  if(dim(ldata)[1] > 0){
    for(j in 1:dim(ldata)[1]){
      plotdata = ldata[j,]
      # Plot a rectangle of the loxal hit
      rect(plotdata$QryStart/scaling,ymin,plotdata$QryEnd/scaling,ymin+(ymax-ymin)*ysplit,col=chromCol(hit)) 
      text(0,(ymax+ymin)/2,hitctg,pos=2,xpd=TRUE,cex=0.8,col=chromCol(hit))
      if(plotdata$SbjStart==1){
        points(plotdata$QryStart/scaling,(ymax+ymin)/2,pch=23,bg=chromCol(hit))
      }
      if(plotdata$SbjEnd==seqData(hit,"Length")){
        points(plotdata$QryEnd/scaling,(ymax+ymin)/2,pch=23,bg=chromCol(hit))
      }
      # Plot the Xdepth data
      if(dim(xdata)[1] > 0 & max(xdata$Pos) >= plotdata$SbjStart){
        xx = xdata[xdata$Pos >= plotdata$SbjStart,]$Pos
        xy = xdata[xdata$Pos >= plotdata$SbjStart,]$X
        if(xx[1] > plotdata$SbjStart){
          xx = c(xx[1],xx)
          xy = c(xy[1],xy)
        }
        xy = c(xy[xx < plotdata$SbjEnd],xy[xx >= plotdata$SbjEnd][1])
        xx = c(xx[xx < plotdata$SbjEnd],xx[xx >= plotdata$SbjEnd][1])
        # Scale xx to Qry positions
        qrange = plotdata$QryEnd-plotdata$QryStart
        srange = plotdata$SbjEnd-plotdata$SbjStart
        qscale = qrange/srange
        (xx = ((xx - plotdata$SbjStart) * qscale) + plotdata$QryStart)
        # Scale xy to yscale
        yheight = ((ymax - ymin) * (1 - ysplit)) / (maxdep * 1.25)
        (xy = (xy * yheight) + ymin + ysplit*(ymax - ymin) )
        # xx and xy should now span from SbjStart to SbjEnd, regardless of whether either is in X Pos data
        lines(xx/scaling,xy,col="blue",type="l",lwd=2)
        (avy = (avdep*yheight) + ymin + ysplit*(ymax - ymin) )
        abline(h=avy,col="red")
        (maxy = (maxdep*yheight) + ymin + ysplit*(ymax - ymin) )
        abline(h=maxy,col="black")
      }else{
        text((plotdata$QryStart+plotdata$QryEnd)/(2*scaling),ymax,"No depth data",pos=1)
      }
    }
  }
}



# This plots a set of contigs or chromosomes against a query using local data and ctype depth data
locAlnDepthPlot = function(qry,hitlist,locdata,ctype="Assembly"){
  # Setup Empty Plot Region and Labels. chromPlotSetup will scale x-axis automatically.
  ymax = length(hitlist)
  ylab = paste(ctype, "hits and read depth")
  ylab = ""
  main = paste(qry,"hit stack read depth")
  chromPlotSetup(xmax=seqData(qry,"Length"),ymax=ymax,ylab=ylab,main=main,ygrid=c(1/6,1.0),yaxt="n")

  # Generate contig plots
  depdata = depthdb[[ctype]]
  xcovdata = xcovdb[[ctype]]
  (medianx = median(xcovdata$MedianX))
  
  i = ymax
  for(hit in hitlist){
    i = i - 1
    (maxdep = max(depdata[depdata$Locus == hit,]$X))
    maxdep = max(maxdep,medianx)
    locXPlot(qry,hit,locdata,depdata,maxdep,medianx,i,i+1,ysplit=0.2)
  }
}

