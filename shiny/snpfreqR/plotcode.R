#i# See the `main.R` file for version, contact and license information.
#i# The `main.R` script loads libraries and contains the initial parameter settings and functions.
#i# This script is loaded by `main.R` and functions as a convenient place for the primary plotting function.
#i# This in primarily to make it easier to find and edit the plotting code as new functions are added.

oldstyle = TRUE

## ~ Plot changes in SNP frequencies across chromosome ~~~~~~~~~~~~~~~~ ##
#># snpdata has the SNP Frequency data
#># chromdata has the chromosome size etc.
#># colbychrom=TRUE uses the colour of the matching chromosome for the points and lines
#># colbychrom=FALSE uses the red for up and blue for down. (Purple for no change.)
chromSNPFreqPlot = function(appdata,chromdata,chrom,alt="Alt",comparison=settings$basename,ymin=0,ymax=1.0,xmin=0,xmax=-1,colbychrom=TRUE,region=FALSE, plotcol="Direction"){
  #DEBUG# writeLines("chromSNPFreqPlot")
  #!# NOTE: This needs tidying a bit now that it is divorced from the old SNPFreqR code.
  
  ### >>> SELECT DATA TO PLOT >>> ###
  #DEBUG# writeLines(appdata$freqlist)
  if(appdata$freqlist == "Loaded data"){
    snpdata = appdata$freqtable
  }else{
    snpdata = appdata$freqdb[[appdata$freqlist]]
  }
  #i# Special TimeLine data
  if(appdata$plothist %in% c("TimeLine","TimeLine (all chrom)")){
    snpdata = appdata$timedb
  }

  ### >>> FILTER SNPs >>> ###
  fdata = filterSNPs(appdata,snpdata)
  snpdata = fdata$snpdata
  
  #desc = fdata$desc
  desc = c(chrom," - ")
  #>> FILTER Description >>#
  #1# Parents
  if(appdata$parlist != "No Filter"){
    desc = c(desc,"Par:",appdata$parlist)
  }
  if(substr(appdata$parlist,1,3) == 'mbg'){   # Single Parent
    if(appdata$paruniq){  # Unique to parent
      if(appdata$parinv){
        desc = c(desc," (UniqInv)")
      }else{
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in parent
      if(appdata$parinv){
        desc = c(desc," (Inv)")
      }
    }
  }
  if(appdata$parlist != "No Filter"){ desc = c(desc,"; ") }
  #2# MBG non-parent Strains
  if(appdata$strlist != "No Filter"){
    desc = c(desc,"MBG1:",appdata$strlist)
  }
  if(substr(appdata$strlist,1,3) == 'mbg'){   # Single strain
    if(appdata$struniq){  # Unique to strain
      if(appdata$strinv){
        desc = c(desc," (UniqInv)")
      }else{
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in strain
      if(appdata$strinv){
        desc = c(desc," (Inv)")
      }
    }
  }
  if(appdata$strlist != "No Filter"){ desc = c(desc,"; ") }
  if(appdata$str2list != "No Filter"){
    desc = c(desc,"MBG2:",appdata$str2list)
  }
  if(substr(appdata$str2list,1,3) == 'mbg'){   # Single strain
    if(appdata$str2uniq){  # Unique to strain
      if(appdata$str2inv){
        desc = c(desc," (UniqInv)")
      }else{
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in strain
      if(appdata$str2inv){
        desc = c(desc," (Inv)")
      }
    }
  }
  if(appdata$str2list != "No Filter"){ desc = c(desc,"; ") }
  #3# Populations
  if(appdata$poplist != "No Filter"){
    desc = c(desc,"Pop:",appdata$poplist)
  }
  if(substr(appdata$poplist,1,3) == 'Pop'){   # Single population
    if(appdata$popuniq){  # Unique to population
      if(appdata$popinv){
        desc = c(desc," (UniqInv)")
      }else{
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in population
      if(appdata$popinv){
        desc = c(desc," (Inv)")
      }
    }
  }
  if(appdata$poplist != "No Filter"){ desc = c(desc,"; ") }
  if(appdata$findft != ""){
    desc = c(desc," (Features)")
  }
  plotfiltdesc = c("[",appdata$freqfilter)
  plotfiltdesc = c(plotfiltdesc,paste0(c("|","+","-","=","|")[c(TRUE,appdata$plotpos,appdata$plotneg,appdata$plotfix,TRUE)],collapse=""))
  plotfiltdesc = c(plotfiltdesc,"]")
  plotfiltdesc = paste0(plotfiltdesc,collapse="")
  
  #4# SNP Types and Effects  
  pop1 = colnames(snpdata)[6]
  if(pop1=="N.Pop00"){ pop1="Pop00" }
  pop2 = colnames(snpdata)[8]
  desc = paste0(c(desc,"SNP:(",appdata$snplist,"|",appdata$efflist,") ",plotfiltdesc),collapse="")
  if(appdata$plothist %in% c("Scatterplot")){
    desc = paste0(desc," ",pop1," to ",pop2," (n = ",nrow(snpdata[snpdata$Locus==chrom,])," / ",nrow(snpdata),")")
  }
  if(appdata$plothist %in% c("Histogram","Histogram (all chrom)")){
    if(appdata$plothist %in% c("Histogram (all chrom)")){
      desc = paste0(desc," ",appdata$plotfield," (n = ",nrow(snpdata),")")
    }else{
      desc = paste0(desc," ",appdata$plotfield," (n = ",nrow(snpdata[snpdata$Locus==chrom,])," / ",nrow(snpdata),")")
    }
  }
  if(appdata$plothist %in% c("TimeLine","TimeLine (all chrom)")){
    if(appdata$plothist %in% c("TimeLine (all chrom)")){
      desc = paste0(desc," ",appdata$plotfield," (n = ",nrow(snpdata),")")
    }else{
      desc = paste0(desc," ",appdata$plotfield," (n = ",nrow(snpdata[snpdata$Locus==chrom,])," / ",nrow(snpdata),")")
    }
  }
  writeLines(desc)
  appdata$intro = desc
  #!# Update to use this: desc = fdata$desc
  #writeLines(desc)
  #appdata$intro = desc
  
  chromdata = appdata$feattable
  comparison = paste0(appdata$complabel," ",pop1," to ",pop2," (n = ",nrow(snpdata[snpdata$Locus==chrom,])," / ",nrow(snpdata),")")
  xmin = appdata$xmin
  xmax = appdata$xmax
  ymin = appdata$ymin
  ymax = appdata$ymax
  
  ### ~ SETUP PLOT ~ ##
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  if(appdata$plothist %in% c("Histogram (all chrom)","TimeLine (all chrom)")){
    plotdata = snpdata
  }else{
    if(alt == "Alt"){
      (plotdata = snpdata[snpdata$AltLocus==chrom,])
      plotdata = plotdata[(order(plotdata$AltPos)),]
    }else{
      (plotdata = snpdata[snpdata$Locus==chrom,])
      plotdata = plotdata[(order(plotdata$Pos)),]
    }
  }
  if(xmax < xmin){
    xmax = max(chromdata[chromdata$Chrom==chrom,]$End)/scaling
  }
  # Setup Title
  if(region == FALSE){ region = chrom }
  if(appdata$titlestyle == "Standard"){
    ptitle = paste(region,comparison)
  }
  if(appdata$titlestyle == "Pure"){
    ptitle = appdata$complabel
  }
  if(appdata$titlestyle == "Info"){
    ptitle = desc
  }
  if(appdata$titlestyle == "None"){
    ptitle = ""
  }
  

  # Setup Plot
  if(appdata$plotfield == "SNPFreq Change"){
    if(ymin == 0){ ymin = -1 }
  }
  
  writeLines(appdata$plothist)
  if(appdata$plothist %in% c("Histogram","Histogram (all chrom)")){
    if(dim(plotdata)[1] < 1){ return(paste0(desc,": No data to plot!")) }
    #Histogram test code = good. Tidy up and add!
    # rgb(1,0,0) = red
    # rgb(0,1,0) = green
    # rgb(0,0,1) = blue
    # rgb(1,1,0) = yellow
    pcol = c()
    for(i in c(0:99)/99){
      pcol = c(pcol,rgb(0,1-i,i))
    }
    pcol = c(pcol,rgb(0,0,0))
    for(i in c(0:99)/99){
      pcol = c(pcol,rgb(1,i,0))
    }
    if(appdata$plotfield %in% c("SNPFreq","Ending SNPFreq")){
      hist(plotdata$MajFreq,breaks=c(-100.5:100.5)/100,xlim=c(ymin,ymax),xlab=appdata$plotfield,col=pcol,main=ptitle)
    }
    if(appdata$plotfield == "Starting SNPFreq"){
      hist(plotdata$MajFreq-plotdata$MajDiff,breaks=c(-100.5:100.5)/100,xlim=c(ymin,ymax),xlab=appdata$plotfield,col=pcol,main=ptitle)
    }
    if(appdata$plotfield == "SNPFreq Change"){
      hist(plotdata$MajDiff,breaks=c(-100.5:100.5)/100,xlim=c(ymin,ymax),xlab=appdata$plotfield,col=pcol,main=ptitle)
    }
    if(appdata$plotfield == "Absolute SNPFreq Change"){
      hist(abs(plotdata$MajDiff),breaks=c(-100.5:100.5)/100,xlim=c(ymin,ymax),xlab=appdata$plotfield,col=pcol,main=ptitle)
    }
    #
    return(desc)
  }
  
  # Set up scatterplot
  if(appdata$plothist %in% c("Scatterplot")){
    plot(c(xmin,xmax),c(ymin,ymax),col="red",type="n",ylab=appdata$plotfield,main=ptitle,xlab=paste("Chromosome Position",units),xaxs="i")
    plotGridlines(xmin,xmax,ymin,ymax,10000/scaling,50000/scaling,0.1,0.5)
  }
  # Set up TimeLine plot
  if(appdata$plothist %in% c("TimeLine","TimeLine (all chrom)")){
    plot(c(0,13),c(ymin,ymax),col="red",type="n",ylab=appdata$plotfield,main=ptitle,xlab=paste("Time Point (Population)"),xaxs="i")
    plotGridlines(0,13,ymin,ymax,1,2,0.1,0.5)
  }

  
  writeLines("Base Plotted")
  
  #else{ plotcol = "Direction" }
  #return(1)
  if(dim(plotdata)[1] < 1){ return(desc) }
  
  
  
  #### >>> DEFINE COLOURS AS NEW FIELDS >>> ###
  plotcol = appdata$plotcol
  if(colbychrom){ plotcol = "Chromosome" }
  if(! oldstyle){
    #i# Default colours are semi-transparent grey, and no fill:
    plotdata$pcol = rgb(0.8,0.8,0.8,0.5)
    plotdata$bg = NA
    plotdata$plwd = 1
    plotdata$ppch = 23
    #i# Colour by parent, or strain
    if(plotcol == "Parent"){
      #transform(plotdata,pcol=ifelse(Parents %in% allpar,appdata$col[[Parents]],pcol))
      for(pparent in allpar){
        if(sum(as.character(plotdata$Parents)==pparent) > 0){
          plotdata[as.character(plotdata$Parents)==pparent,]$pcol = appdata$col[[pparent]]
        }
      }
    }
    if(plotcol == "MBG Strain"){
      for(pstrain in allmbg){
        plotdata[as.character(plotdata$Strains)==pstrain,]$pcol = appdata$col[[pstrain]]
      }
    }
    if(plotcol == "Chromosome"){
      for(chrom in appdata$chromlist){
        plotdata[as.character(plotdata$Locus)==chrom,]$pcol = appdata$col[[chrom]]
      }
    }
    #i# Colour by direction of change
    if(! plotcol %in% c("Parent","MBG Strain","Chromosome")){
      
      changedata = plotdata$MajDiff > -2 # All TRUE
      #i# MBG Strain 1
      if(plotcol == "Strain 1 Selection"){ changedata = regexpr(appdata$strlist,plotdata$Strains) >= 0 }
      # Colour mutations only
      if(plotcol == "No Parent"){ changedata = as.character(plotdata$Parents) == "" }
      # Colour fixation only
      if(plotcol == "Fixation"){ changedata = snpdata$MajFreq >= (1.0-appdata$fixfreq) | plotdata$MajFreq <= appdata$fixfreq }
      # SNPType only
      if(plotcol == "SNPType Selection"){
        if(appdata$snplist == "CDS"){ changedata = plotdata$SNPType != "SNP" }
        else{ changedata = plotdata$SNPType == appdata$snplist }
      }
      # SNPEffect only
      if(plotcol == "SNPEffect Selection"){ changedata = plotdata$SNPEffect == appdata$efflist }

      if(sum(plotdata$MajDiff == 0 & changedata) > 0){
        plotdata[plotdata$MajDiff == 0 & changedata,]$pcol = rgb(1,0.5,1)
      }
      for(r in which(plotdata$MajDiff > 0 & changedata)){
        i = min(1,max(0.2,2*plotdata[r,]$MajDiff))
        plotdata[r,]$pcol = rgb(1,0,0,i)
      }
      for(r in which(plotdata$MajDiff < 0 & changedata)){
        i = max(-1.0,min(-0.2,2*plotdata[r,]$MajDiff))
        plotdata[r,]$pcol = rgb(0,0,1,-i)
      }
    }

    writeLines(paste(plotcol,"Colours set"))
  }

  ### >>> Set plotting shape >>> ###
  plotdata$ppch = 23
  if(sum(plotdata$MajDiff > 0)){
    plotdata[plotdata$MajDiff > 0,]$ppch = 24
  }
  if(sum(plotdata$MajDiff < 0)){
    plotdata[plotdata$MajDiff < 0,]$ppch = 25
  }
  
  ### >>> REORDER DATA FOR PLOTTING OVERLAYS >>> ###
  if(appdata$plothist %in% c("TimeLine","TimeLine (all chrom)")){
    #snpdata = appdata$timedb
    #Try ordering. Might make this an option
    #snpdata = snpdata[order(snpdata$Pop00),]
    #snpdata = snpdata[order(snpdata$Pop13),]
    
    #plotdata = plotdata[order(plotdata$Pop00),]
    #plotdata = plotdata[order(plotdata$Pop13),]
    
    #!# Rather than colours, create a sorting field
    
    #!# Reorder based on colouring! (Grey first, then rest)
    if(! oldstyle){ #!# Replace with appdata$oldstyle
      plotdata
    }
  }
  
  
  
  # Plot Data
  plotfield = appdata$plotfield
  if(appdata$plothist %in% c("TimeLine","TimeLine (all chrom)")){ plotfield = "TimeLine" }
  for(xi in 1:length(plotdata$Pos)){
    pdata = plotdata[xi,]
    
    #Check whether to plot
    #if(pdata$MajDiff > 0){
    #  if(appdata$plotpos == FALSE){ next }
    #}else if(pdata$MajDiff == 0){
    #  if(appdata$plotfix == FALSE){ next }
    #}else{
    #  if(appdata$plotneg == FALSE){ next }
    #}
    
    # Set plot symbol
    if(pdata$MajDiff > 0){
      ppch=24
    }else if(pdata$MajDiff == 0){
      ppch=23
    }else{
      ppch=25
    }
    
    ppch = pdata$ppch
    
    if(oldstyle){ #Use old style colours if true
      ### Set plot colours
      ## Colour line and outlines
      # Colour by chromosome
      if(plotcol == "Chromosome"){
        if(alt == "Alt"){
          pchrom = as.character(pdata$Locus)
          pcol = settings$col[[pchrom]]
        }else{
          pchrom = as.character(pdata$AltLocus)
          pcol = settings$col[[pchrom]]
        }
      }
      # Colour by parent
      pparent = as.character(pdata$Parents)
      if(plotcol == "Parent"){
        if(pparent %in% names(appdata$col)){
          pcol = appdata$col[[pparent]]
        }else{
          pcol = rgb(0.8,0.8,0.8,0.5)
        }
      }
      # Colour by MBG Strain
      pstrain = as.character(pdata$Strains)
      if(plotcol == "MBG Strain"){
        if(pstrain %in% names(appdata$col)){
          pcol = appdata$col[[pstrain]]
        }else{
          pcol = rgb(0.8,0.8,0.8,0.5)
        }
      }
      # Colour by direction of change
      if(plotcol %in% c("Direction","Strain 1 Selection","No Parent","Fixation","SNPType Selection","SNPEffect Selection")){
        if(pdata$MajDiff > 0){
          i = min(1,max(0.2,2*pdata$MajDiff))
          pcol = rgb(1,1-i,1-i)
          pcol = rgb(1,0,0,i)
        }else if(pdata$MajDiff == 0){
          pcol = rgb(1,0.5,1)
          pcol = rgb(1,0.5,1,0.5)
        }else{
          i = max(-1.0,min(-0.2,2*pdata$MajDiff))
          pcol = rgb(1+i,1+i,1)
          pcol = rgb(0,0,1,-i)
        }
      }
      
      # Colour by direction of change
      if(plotcol == "Strain 1 Selection"){
        if(regexpr(appdata$strlist,pdata$Strains) < 0){
          pcol = rgb(0.8,0.8,0.8,0.5)
        }
      }
      # Colour mutations only
      if(plotcol == "No Parent"){
        if(pparent != ""){
          pcol = rgb(0.8,0.8,0.8,0.5)
        }
      }
      # Colour fixation only
      if(plotcol == "Fixation"){
        if(pdata$MajFreq < (1.0-appdata$fixfreq) & pdata$MajFreq > appdata$fixfreq){
          pcol = rgb(0.8,0.8,0.8,0.5)
        }
      }
      # SNPType only
      if(plotcol == "SNPType Selection"){
        if(appdata$snplist == "CDS"){
          if(pdata$SNPType == "SNP"){  
            pcol = rgb(0.8,0.8,0.8,0.5)
          }
        }else{
          if(pdata$SNPType != appdata$snplist){
            pcol = rgb(0.8,0.8,0.8,0.5)
          }
        }
      }
      # SNPEffect only
      if(plotcol == "SNPEffect Selection"){
        if(pdata$SNPEffect != appdata$efflist){
          pcol = rgb(0.8,0.8,0.8,0.5)
        }
      }
    }else{
      pcol = pdata$pcol
    }
    
    ## Background symbol highlights
    pbg = NA
    plwd = 1
    if(pdata$MajProb <= appdata$pcut & appdata$plothightlights & pcol != rgb(0.8,0.8,0.8,0.5)){ #    plotcol %in% c("Direction","Strain 1 Selection","No Parent")){
      p = pdata$MajProb
      i = max(0,(10 + log10(p))/10.0)
      if(pdata$MajDiff >= 0){
        pbg = rgb(1,1,i)
      }else{
        pbg = rgb(i,1,i)
      }  
      #plwd = 2 - i
    }
    pparent = as.character(pdata$Parents)
    if(appdata$plothightlights & plotcol == "Parent" & pparent %in% names(appdata$col)){ 
      pbg = pcol 
    }
    pstrain = as.character(pdata$Strains)
    if(appdata$plothightlights & plotcol == "MBG Strain" & pstrain %in% names(appdata$col)){ 
      pbg = pcol 
    }
    
    
    # TimeLine plot
    if(plotfield == "TimeLine"){
      # Plot line of frequencies
      xdata = c(0,4,6,9,10,11,12,13)
      ydata = c(pdata$Pop00,pdata$Pop04,pdata$Pop06,pdata$Pop09,pdata$Pop10,pdata$Pop11,pdata$Pop12,pdata$Pop13)
      lines(xdata,ydata,col=pcol,type="l",lwd=plwd)
    }
    # "SNPFreq" = Basic plot
    if(plotfield == "SNPFreq"){
      # Plot lines and symbols
      if(alt == "Alt"){
        lines(c(pdata$AltPos/scaling,pdata$AltPos/scaling),c(pdata$MajFreq,pdata$MajFreq-pdata$MajDiff),col=pcol,type="l",lwd=plwd)
        points(c(pdata$AltPos/scaling),c(pdata$MajFreq),col=pcol,pch=ppch,bg=pbg)
      }else{
        lines(c(pdata$Pos/scaling,pdata$Pos/scaling),c(pdata$MajFreq,pdata$MajFreq-pdata$MajDiff),col=pcol,type="l",lwd=plwd)
        points(c(pdata$Pos/scaling),c(pdata$MajFreq),col=pcol,pch=ppch,bg=pbg)
      }
    }
    if(plotfield == "Starting SNPFreq"){
      # Plot lines and symbols
      if(alt == "Alt"){
        points(c(pdata$AltPos/scaling),c(pdata$MajFreq-pdata$MajDiff),col=pcol,pch=ppch,bg=pbg)
      }else{
        points(c(pdata$Pos/scaling),c(pdata$MajFreq-pdata$MajDiff),col=pcol,pch=ppch,bg=pbg)
      }
    }
    if(plotfield == "Ending SNPFreq"){
      # Plot lines and symbols
      if(alt == "Alt"){
        points(c(pdata$AltPos/scaling),c(pdata$MajFreq),col=pcol,pch=ppch,bg=pbg)
      }else{
        points(c(pdata$Pos/scaling),c(pdata$MajFreq),col=pcol,pch=ppch,bg=pbg)
      }
    }
    #"SNPFreq Change" = centre P1 freq to be 0.0
    if(plotfield == "SNPFreq Change"){
      if(ymin == 0){ ymin = -1 }
      # Plot lines and symbols
      if(alt == "Alt"){
        lines(c(pdata$AltPos/scaling,pdata$AltPos/scaling),c(0,pdata$MajDiff),col=pcol,type="l",lwd=plwd)
        points(c(pdata$AltPos/scaling),c(pdata$MajDiff),col=pcol,pch=ppch,bg=pbg)
      }else{
        lines(c(pdata$Pos/scaling,pdata$Pos/scaling),c(0,pdata$MajDiff),col=pcol,type="l",lwd=plwd)
        points(c(pdata$Pos/scaling),c(pdata$MajDiff),col=pcol,pch=ppch,bg=pbg)
      }
    }
    #"Absolute SNPFreq Change" = centre P1 freq to be 0.0 and invert negative changes
    if(plotfield == "Absolute SNPFreq Change"){
      # Plot lines and symbols
      if(alt == "Alt"){
        lines(c(pdata$AltPos/scaling,pdata$AltPos/scaling),c(0,abs(pdata$MajDiff)),col=pcol,type="l",lwd=plwd)
        points(c(pdata$AltPos/scaling),c(abs(pdata$MajDiff)),col=pcol,pch=ppch,bg=pbg)
      }else{
        lines(c(pdata$Pos/scaling,pdata$Pos/scaling),c(0,abs(pdata$MajDiff)),col=pcol,type="l",lwd=plwd)
        points(c(pdata$Pos/scaling),c(abs(pdata$MajDiff)),col=pcol,pch=ppch,bg=pbg)
      }
    }
  }#lines(plotdata$Pos/scaling,plotdata$ContigNum-plotdata$ChromNum,col="blue",type="l")
  writeLines("Plot finished.")
  return(desc)
}
