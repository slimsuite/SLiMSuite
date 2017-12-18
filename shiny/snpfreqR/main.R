################# ::: APP INFO ::: ######################
info = list(
  apptitle = "SNPFreqR MBG Evolution",
  version = "0.7.0",
  lastedit = "12/12/2017",
  author = "Richard J. Edwards",
  contact = "richard.edwards@unsw.edu.au",
  description = "Custom SNP Frequency analysis for MBGEvolution project."
)
#i# The main.R script load libraries and contains the initial parameter settings and functions.
#i# This script is used by the `ui.R` and `server.R` files.

################# ::: HISTORY ::: ######################
# V0.1.0 - Development version.
# V0.2.0 - Basic functions with filtering.
# V0.3.0 - Added SNP count outputs.
# V0.4.0 - Added SNPMap table and extra title options.
# V0.4.1 - Bug fixes with Pop filter and extra N stats.
# V0.5.0 - Added multiple SNPFreq files and improved PNG output.
# V0.6.0 - Added locus_tag filer.
# V0.6.1 - Fixed PNG file name bug and moved other RJE R scripts to main directory for ease of sharing.
# V0.7.0 - Updated interface for Microbiogen, plus initial in-app documentation.

################## ::: LICENSE ::: #####################
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# SLiMEnrich program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License (gpl.md)
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

################# ::: KNOWN BUGS ::: ######################
#!# List bugs here.
#!# Might need to run source("main.R")

################# ::: TO DO LIST ::: ######################
# [ ] : Filter the Features table.


############### ::: GENERAL SETUP ::: ##################
#i# This used to be the setup.R code...

############### ::: REQUIRED LIBRARIES ::: ##################
# Check whether packages of interest are installed
is_installed = function(mypkg) is.element(mypkg, installed.packages()[,1]) 
# Install library if not already installed
# Run a for-loop of all the package names listed below in the function call
# with the list of packages: load_or_install(c("pkg1", "pkg2",..., "pkgn"))
load_or_install = function(package_names) 
{ 
  for(package_name in package_names) 
  { 
    if(!is_installed(package_name)) 
    { 
      #install.packages(package_name,repos="http://lib.stat.cmu.edu/R/CRAN") 
      install.packages(package_name)
    } 
    library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE) 
  } 
}
# Load or install required libraries
load_or_install(c("shiny", "httr", "DT", "markdown","plyr", "shinyFiles"))

############### ::: LOAD R LIBRARIES ::: ##################
rdir = './' # Was: '../../libraries/r/'
source(paste(rdir,"rje_col.r",sep=""))
source(paste(rdir,"rje_misc.r",sep=""))
source(paste(rdir,"rje_genomics.r",sep=""))


############### ::: SET DEFAULTS ::: ##################
settings = list(
  # Setup scaling and units
  scaling = 1,
  units = "kb",
  features = TRUE,
  testing=FALSE,
  fdrplot = FALSE,
  
  inputpath = "../../analysis/2017-11-27.SNPFreq",
  #inputpath = "/Users/redwards/OneDrive - UNSW/projects/MBGEvolution-Feb17/analysis/2017-11-27.SNPFreq",
  pngpath = "../../analysis/2017-11-27.SNPFreq/plots",
  
  
  #i# URL of REST servers
  resturl = "http://rest.slimsuite.unsw.edu.au/",
  #i# URL of Dev REST servers
  devurl = "http://restdev.slimsuite.unsw.edu.au/",
  #i# Standard REST server &rest keys
  stdkeys = c("status", "version", "ini", "log", "warnings", "errors", "outfmt", "help", "jobid", "intro", "prog"),
  #i# Subset of REST server keys that are not returned by restKeys()
  restkeys = c("errors", "outfmt", "help", "jobid", "intro", "prog"),
  #i# Define REST output formats
  csv = c("main","occ"),
  tsv = c(),
  #i# Whether to run in debugging mode
  debug = FALSE,
  #i# Default program
  prog = "zen",
  #i# Default JobID
  jobid = "17092300003"   # Test with 17092500002
)

if(getwd() == "/Users/redwards/Dropbox/_Repository_/slimsuite/shiny/snpfreqR"){
  settings$inputpath = "/Users/redwards/OneDrive - UNSW/projects/MBGEvolution-Feb17/analysis/2017-11-27.SNPFreq"
  settings$pngpath = "/Users/redwards/OneDrive - UNSW/projects/MBGEvolution-Feb17/analysis/2017-11-27.SNPFreq/plots"
}

## ~ Set up global plot settings (SNPFreq) ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
settings$units = "kb"
settings$pngwidth = 2400
settings$pngheight = 1600
settings$dualpng = 1600
settings$pointsize = 36
settings$pcut = 0.01
settings$fdrplot = FALSE
settings$multiplot = FALSE

## ~ Setup Chromosome plot settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
settings$chrom = NA
settings$chrom = "sgdI_YEAST__BK006935" # "sgdIII_YEAST__BK006937"
settings$ymin=0
settings$ymax=1.0
settings$xmin=0
settings$xmax=-1
settings$plotpos = TRUE
settings$plotneg = TRUE
settings$plotheight = 750

# Update 
if(settings$units == "kb"){
  settings$scaling = 1000
}
if(settings$units == "Mb"){
  settings$scaling = 1e6
}
# Setup colour profile
settings = featureDefaults(settings)
#if(settings$features){
#  writeLines("Features for plotting:")
#  writeLines(settings$plotft)
#}
settings = featureColours(settings)


############### ::: SETUP DATA ::: ##################
#i# This function is called by server.R in a reactiveValues() call.
#i# This way, it will be updated whenever the settings that affect it are changed.
setupData = function(){
  emptydb = list()
  for(rkey in c("freqtable","feattable")){
    emptydb[[rkey]] = data.frame()
  }
  emptydb$freqdb = list()  # List for multiple SNP data.frames
  emptydb$freqfiles = c()  # Vector of filenames for snpdb
  emptydb$freqlist = ""    # Selected population comparison to plot.
  emptydb$freqloaded = FALSE   # Stores whether SNP frequencies have been loaded.
  emptydb$featloaded = FALSE   # Stores whether features have been loaded.
  emptydb$loads = 0      # This counter is to prevent action buttons from self-triggering
  emptydb$pngsaves = 0   # This counter is to prevent action buttons from self-triggering
  emptydb$intro = "Load SNP Frequencies and Reference Features. (See Instructions Tab)"
  emptydb$settings = settings
  # Return data list
  return(emptydb)
}



#i#featurePlot("Reference","BK006934",xmin=1,xmax=300,acconly=FALSE)
#i# NOTE: featurePlot() uses pre-scaled xmin and xmax. 

############# ::: DEFINE SNPFREQ METHODS ::: ##################

## ~ FreqTable formatting cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
freqTableCleanup = function(freqdb){
  if(! 'AltPos' %in% colnames(freqdb)){
    freqdb$AltPos = freqdb$Pos
    freqdb$AltLocus = freqdb$Locus
  }
  freqdb$Pops = as.character(freqdb$Pops)
  if(length(freqdb[freqdb$Pops == "Pops",]$Pops) > 0){
    freqdb[freqdb$Pops == "Pops",]$Pops = ""
  }
  freqdb$Pops = as.factor(freqdb$Pops)
  freqdb$SNPEffect = as.character(freqdb$SNPEffect)
  if(length(freqdb[freqdb$SNPType == "NS",]$SNPEffect) > 0){
    freqdb[freqdb$SNPType == "NS",]$SNPEffect = "Nonsyn"
  }
  if(length(freqdb[freqdb$SNPType == "SYN",]$SNPEffect) > 0){
    freqdb[freqdb$SNPType == "SYN",]$SNPEffect = "Syn"
  }
  if(length(freqdb[freqdb$SNPType == "NON",]$SNPEffect) > 0){
    freqdb[freqdb$SNPType == "NON",]$SNPEffect = "STOP"
  }
  if(length(freqdb[freqdb$SNPType == "EXT",]$SNPEffect) > 0){
    freqdb[freqdb$SNPType == "EXT",]$SNPEffect = "Extension"
  }
  freqdb$SNPEffect = as.factor(freqdb$SNPEffect)
  return(freqdb)
}


## ~ FreqTable SNP filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
filterSNPs = function(appdata,freqdb){
  snpdata = freqdb
  ### Subsetting data
  if(appdata$freqfilter == "Fixed"){
    snpdata = snpdata[snpdata$MajFreq >= (1.0-appdata$fixfreq),]  
  }
  if(appdata$freqfilter == "Lost"){
    snpdata = snpdata[snpdata$MajFreq <= appdata$fixfreq,]  
  }
  if(appdata$freqfilter == "Polymorphic"){
    snpdata = snpdata[snpdata$MajFreq >= appdata$fixfreq & snpdata$MajFreq <= (1.0-appdata$fixfreq),]  
  }
  desc = c(appdata$freqfilter," ")
  #1# Parents
  desc = c(desc,"Par:",appdata$parlist)
  if(substr(appdata$parlist,1,3) == 'mbg'){   # Single Parent
    if(appdata$paruniq){  # Unique to parent
      if(appdata$parinv){
        snpdata = snpdata[snpdata$Parents != appdata$parlist,]  
        desc = c(desc," (UniqInv)")
      }else{
        snpdata = snpdata[snpdata$Parents == appdata$parlist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in parent
      if(appdata$parinv){
        snpdata = snpdata[regexpr(appdata$parlist,snpdata$Parents) < 0,]  
        desc = c(desc," (Inv)")
      }else{
        snpdata = snpdata[regexpr(appdata$parlist,snpdata$Parents) > 0,]  
      }
    }
  }
  if(appdata$parlist == "Any"){   # Any parent
    if(appdata$parinv){
      snpdata = snpdata[regexpr("mbg",snpdata$Parents) < 0,]  
    }else{
      snpdata = snpdata[regexpr("mbg",snpdata$Parents) > 0,]  
    }
  }  
  if(appdata$parlist == "None"){  # No parent
    snpdata = snpdata[snpdata$Parents == "",]  
  }
  
  
  #2# MBG non-parent Strains
  desc = c(desc,"; MBG:",appdata$strlist)
  if(substr(appdata$strlist,1,3) == 'mbg'){   # Single strain
    if(appdata$struniq){  # Unique to strain
      if(appdata$strinv){
        snpdata = snpdata[snpdata$Strains != appdata$strlist,]  
        desc = c(desc," (UniqInv)")
      }else{
        snpdata = snpdata[snpdata$Strains == appdata$strlist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in strain
      if(appdata$strinv){
        snpdata = snpdata[regexpr(appdata$strlist,snpdata$Strains) < 0,]  
        desc = c(desc," (nv)")
      }else{
        snpdata = snpdata[regexpr(appdata$strlist,snpdata$Strains) > 0,]  
      }
    }
  }
  if(appdata$strlist == "Any"){   # Any strain
    if(appdata$strinv){
      snpdata = snpdata[regexpr("mbg",snpdata$Strains) < 0,]  
    }else{
      snpdata = snpdata[regexpr("mbg",snpdata$Strains) > 0,]  
    }
  }  
  if(appdata$strlist == "None"){  # No strain
    snpdata = snpdata[snpdata$Strains == "",]  
  }
  
  #3# Populations
  desc = c(desc,"; Pop:",appdata$poplist)
  if(substr(appdata$poplist,1,3) == 'Pop'){   # Single population
    if(appdata$popuniq){  # Unique to population
      if(appdata$popinv){
        desc = c(desc," (UniqInv)")
        snpdata = snpdata[snpdata$Pops != appdata$poplist,]  
      }else{
        snpdata = snpdata[snpdata$Pops == appdata$poplist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in population
      if(appdata$popinv){
        snpdata = snpdata[regexpr(appdata$poplist,snpdata$Pops) < 0,]  
        desc = c(desc," (Inv)")
      }else{
        snpdata = snpdata[regexpr(appdata$poplist,snpdata$Pops) > 0,]  
      }
    }
  }
  if(appdata$poplist == "Any"){   # Any population
    if(appdata$popinv){
      snpdata = snpdata[regexpr("Pop",snpdata$Pops) < 0,]  
    }else{
      snpdata = snpdata[regexpr("Pop",snpdata$Pops) > 0,]  
    }
  }  
  if(appdata$poplist == "None"){  # No population
    snpdata = snpdata[snpdata$Pops == "",]  
  }
  
  #4# SNP Types and Effects  
  if(appdata$snplist != "All"){
    if(appdata$snplist == "CDS"){
      snpdata = snpdata[snpdata$SNPType != "SNP",]  
    }else{
      snpdata = snpdata[snpdata$SNPType == appdata$snplist,]  
    }
  }
  if(appdata$efflist != "All"){
    snpdata = snpdata[snpdata$SNPEffect == appdata$efflist,]  
  }
  pop1 = colnames(snpdata)[6]
  pop2 = colnames(snpdata)[8]
  desc = paste0(c(desc,"; SNP:(",appdata$snplist,"|",appdata$efflist,")"),collapse="")
  desc = paste0(desc," ",pop1," to ",pop2," (n = ",nrow(snpdata),")")
  writeLines(desc)
  
  chromdata = appdata$feattable
  xmin = appdata$xmin
  xmax = appdata$xmax
  ymin = appdata$ymin
  ymax = appdata$ymax
  
  # Locus_tag restriction
  if(appdata$findft != ""){
    writeLines("locus_tag restriction")
    writeLines(appdata$findft)
    writeLines(">>>")
    ftchrom = as.character(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$Chrom[1])
    #writeLines(ftchrom)
    fmin = min(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$start)
    #writeLines(fmin)
    fmax = max(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$end)
    writeLines("Looking for ftsnp")
    if(ftchrom %in% snpdata$Locus){
      snpdata = snpdata[snpdata$Locus == ftchrom & snpdata$Pos >= fmin & snpdata$Pos <= fmax,]  
      if(nrow(snpdata) < 1){
        return(paste(ftchrom,"not found in filtered SNP data."))
      }
    }else{
      return(paste(ftchrom,"not found in filtered SNP data."))
    }
  }
  appdata$intro = desc
  return(snpdata)
}


## ~ Chromosome copy number plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
chromCNVPlot = function(chromdata,chrom,xmin=0,xmax=-1){
  # Plot of the Hits betweeen Chromosomes and Contigs(+) or Chromosomes(-):
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  (plotdata = chromdata[chromdata$Chrom==chrom,])
  plotdata = plotdata[(order(plotdata$Start)),]
  if(xmax < 0){
    xmax = max(plotdata$End)/scaling
  }
  ymin = 0
  ymax = max(plotdata$Copy)
  # Setup Plot
  plot(c(xmin,xmax),c(ymin,ymax),col="red",type="n",ylab="Copy Number",main=paste(chrom,"CNV"),xlab=paste("Chromosome Position",units),xaxs="i")
  plotGridlines(xmin,xmax,ymin,ymax,10000/scaling,50000/scaling,1,10)
  # Plot Data
  for(xi in 1:length(plotdata$Start)){
    recdata = plotdata[xi,]
    rect(recdata$Start/scaling,0.0,recdata$End/scaling,recdata$Copy,col=settings$col[[chrom]])
  }#lines(plotdata$Pos/scaling,plotdata$ContigNum-plotdata$ChromNum,col="blue",type="l")
}

## ~ Plot changes in SNP frequencies across chromosome ~~~~~~~~~~~~~~~~ ##
#># snpdata has the SNP Frequency data
#># chromdata has the chromosome size etc.
#># colbychrom=TRUE uses the colour of the matching chromosome for the points and lines
#># colbychrom=FALSE uses the red for up and blue for down. (Purple for no change.)
chromSNPFreqPlot = function(appdata,chromdata,chrom,alt="Alt",comparison=settings$basename,ymin=0,ymax=1.0,xmin=0,xmax=-1,colbychrom=TRUE,region=FALSE){
  #!# Add switch for time series plot
  timeplot = FALSE
  writeLines("chromSNPFreqPlot")
  # Plot of the Hits betweeen Chromosomes and Contigs(+) or Chromosomes(-):
  if(timeplot){
    snpdata = appdata$snptable
  }else{
    writeLines(appdata$freqlist)
    if(appdata$freqlist == "Loaded data"){
      snpdata = appdata$freqtable
    }else{
      snpdata = appdata$freqdb[[appdata$freqlist]]
    }
  }
  desc = c(chrom," - ")
  
  #>> FILTER START >>#
  ### Subsetting data
  if(appdata$freqfilter == "Fixed"){
    snpdata = snpdata[snpdata$MajFreq >= (1.0-appdata$fixfreq),]  
  }
  if(appdata$freqfilter == "Lost"){
    snpdata = snpdata[snpdata$MajFreq <= appdata$fixfreq,]  
  }
  if(appdata$freqfilter == "Polymorphic"){
    snpdata = snpdata[snpdata$MajFreq >= appdata$fixfreq & snpdata$MajFreq <= (1.0-appdata$fixfreq),]  
  }
  #1# Parents
  desc = c(desc,"Par:",appdata$parlist)
  if(substr(appdata$parlist,1,3) == 'mbg'){   # Single Parent
    if(appdata$paruniq){  # Unique to parent
      if(appdata$parinv){
        snpdata = snpdata[snpdata$Parents != appdata$parlist,]  
        desc = c(desc," (UniqInv)")
      }else{
        snpdata = snpdata[snpdata$Parents == appdata$parlist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in parent
      if(appdata$parinv){
        snpdata = snpdata[regexpr(appdata$parlist,snpdata$Parents) < 0,]  
        desc = c(desc," (Inv)")
      }else{
        snpdata = snpdata[regexpr(appdata$parlist,snpdata$Parents) > 0,]  
      }
    }
  }
  if(appdata$parlist == "Any"){   # Any parent
    if(appdata$parinv){
      snpdata = snpdata[regexpr("mbg",snpdata$Parents) < 0,]  
    }else{
      snpdata = snpdata[regexpr("mbg",snpdata$Parents) > 0,]  
    }
  }  
  if(appdata$parlist == "None"){  # No parent
    snpdata = snpdata[snpdata$Parents == "",]  
  }
  
    
  #2# MBG non-parent Strains
  desc = c(desc,"; MBG:",appdata$strlist)
  if(substr(appdata$strlist,1,3) == 'mbg'){   # Single strain
    if(appdata$struniq){  # Unique to strain
      if(appdata$strinv){
        snpdata = snpdata[snpdata$Strains != appdata$strlist,]  
        desc = c(desc," (UniqInv)")
      }else{
        snpdata = snpdata[snpdata$Strains == appdata$strlist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in strain
      if(appdata$strinv){
        snpdata = snpdata[regexpr(appdata$strlist,snpdata$Strains) < 0,]  
        desc = c(desc," (nv)")
      }else{
        snpdata = snpdata[regexpr(appdata$strlist,snpdata$Strains) > 0,]  
      }
    }
  }
  if(appdata$strlist == "Any"){   # Any strain
    if(appdata$strinv){
      snpdata = snpdata[regexpr("mbg",snpdata$Strains) < 0,]  
    }else{
      snpdata = snpdata[regexpr("mbg",snpdata$Strains) > 0,]  
    }
  }  
  if(appdata$strlist == "None"){  # No strain
    snpdata = snpdata[snpdata$Strains == "",]  
  }
  
  #3# Populations
  desc = c(desc,"; Pop:",appdata$poplist)
  if(substr(appdata$poplist,1,3) == 'Pop'){   # Single population
    if(appdata$popuniq){  # Unique to population
      if(appdata$popinv){
        desc = c(desc," (UniqInv)")
        snpdata = snpdata[snpdata$Pops != appdata$poplist,]  
      }else{
        snpdata = snpdata[snpdata$Pops == appdata$poplist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in population
      if(appdata$popinv){
        snpdata = snpdata[regexpr(appdata$poplist,snpdata$Pops) < 0,]  
        desc = c(desc," (Inv)")
      }else{
        snpdata = snpdata[regexpr(appdata$poplist,snpdata$Pops) > 0,]  
      }
    }
  }
  if(appdata$poplist == "Any"){   # Any population
    if(appdata$popinv){
      snpdata = snpdata[regexpr("Pop",snpdata$Pops) < 0,]  
    }else{
      snpdata = snpdata[regexpr("Pop",snpdata$Pops) > 0,]  
    }
  }  
  if(appdata$poplist == "None"){  # No population
    snpdata = snpdata[snpdata$Pops == "",]  
  }
  
  #4# SNP Types and Effects  
  if(appdata$snplist != "All"){
    if(appdata$snplist == "CDS"){
      snpdata = snpdata[snpdata$SNPType != "SNP",]  
    }else{
      snpdata = snpdata[snpdata$SNPType == appdata$snplist,]  
    }
  }
  if(appdata$efflist != "All"){
    snpdata = snpdata[snpdata$SNPEffect == appdata$efflist,]  
  }
  pop1 = colnames(snpdata)[6]
  pop2 = colnames(snpdata)[8]
  desc = paste0(c(desc,"; SNP:(",appdata$snplist,"|",appdata$efflist,")"),collapse="")
  desc = paste0(desc," ",pop1," to ",pop2," (n = ",nrow(snpdata[snpdata$Locus==chrom,])," / ",nrow(snpdata),")")
  writeLines(desc)
  
  chromdata = appdata$feattable
  comparison = paste0(appdata$complabel," ",pop1," to ",pop2," (n = ",nrow(snpdata[snpdata$Locus==chrom,])," / ",nrow(snpdata),")")
  xmin = appdata$xmin
  xmax = appdata$xmax
  ymin = appdata$ymin
  ymax = appdata$ymax
  
  # Locus_tag restriction
  if(appdata$findft != ""){
    writeLines("locus_tag restriction")
    writeLines(appdata$findft)
    writeLines(">>>")
    ftchrom = as.character(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$Chrom[1])
    #writeLines(ftchrom)
    fmin = min(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$Start)
    #writeLines(fmin)
    fmax = max(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$End)
    writeLines("Looking for ftsnp")
    if(ftchrom %in% snpdata$Locus){
      snpdata = snpdata[snpdata$Locus == ftchrom & snpdata$Pos >= fmin & snpdata$Pos <= fmax,]  
      if(nrow(snpdata) < 1){
        return(paste(ftchrom,"not found in filtered SNP data."))
      }
    }else{
      return(paste(ftchrom,"not found in filtered SNP data."))
    }
  }
  #<< FILTER END <<#

  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  if(alt == "Alt"){
    (plotdata = snpdata[snpdata$AltLocus==chrom,])
    plotdata = plotdata[(order(plotdata$AltPos)),]
  }else{
    (plotdata = snpdata[snpdata$Locus==chrom,])
    plotdata = plotdata[(order(plotdata$Pos)),]
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
  
  # Time series plot option
  if(timeplot){
    plotdata
    xmin
    xmax
    ymin
    ymax
    ptitle
  }
  
  # Setup Plot
  plot(c(xmin,xmax),c(ymin,ymax),col="red",type="n",ylab="SNP Freq",main=ptitle,xlab=paste("Chromosome Position",units),xaxs="i")
  plotGridlines(xmin,xmax,ymin,ymax,10000/scaling,50000/scaling,0.1,0.5)
  writeLines("Base Plotted")
  #return(1)
  if(dim(plotdata)[1] < 1){ return(0) }
  # Plot Data
  for(xi in 1:length(plotdata$Pos)){
    pdata = plotdata[xi,]
    # Set plot symbol
    if(pdata$MajDiff > 0){
      ppch=24
    }else if(pdata$MajDiff == 0){
      ppch=23
    }else{
      ppch=25
    }
    # Set plot colours
    if(pdata$MajProb <= settings$pcut){
      p = pdata$MajProb
      i = max(0,(10 + log10(p))/10.0)
      pbg = rgb(1,1,i)
      plwd = 2 - i
    }else{
      pbg = NA
      plwd = 1
    }
    if(colbychrom){
      if(alt == "Alt"){
        pchrom = as.character(pdata$Locus)
        pcol = settings$col[[pchrom]]
      }else{
        pchrom = as.character(pdata$AltLocus)
        pcol = settings$col[[pchrom]]
      }
    }else{
      if(pdata$MajDiff > 0){
        if(appdata$plotpos == FALSE){ next }
        i = min(1,max(0.2,2*pdata$MajDiff))
        pcol = rgb(1,1-i,1-i)
      }else if(pdata$MajDiff == 0){
        pcol = rgb(1,0.5,1)
        if(appdata$plotfix == FALSE){ next }
      }else{
        if(appdata$plotneg == FALSE){ next }
        i = max(-1.0,min(-0.2,2*pdata$MajDiff))
        pcol = rgb(1+i,1+i,1)
      }
    }
    # Plot lines and symbols
    if(alt == "Alt"){
      lines(c(pdata$AltPos/scaling,pdata$AltPos/scaling),c(pdata$MajFreq,pdata$MajFreq-pdata$MajDiff),col=pcol,type="l",lwd=plwd)
      points(c(pdata$AltPos/scaling),c(pdata$MajFreq),col=pcol,pch=ppch,bg=pbg)
    }else{
      lines(c(pdata$Pos/scaling,pdata$Pos/scaling),c(pdata$MajFreq,pdata$MajFreq-pdata$MajDiff),col=pcol,type="l",lwd=plwd)
      points(c(pdata$Pos/scaling),c(pdata$MajFreq),col=pcol,pch=ppch,bg=pbg)
    }
  }#lines(plotdata$Pos/scaling,plotdata$ContigNum-plotdata$ChromNum,col="blue",type="l")
  writeLines("Plot finished.")
  return(desc)
}

## ~ Plot changes in SNP frequencies across chromosome, assuming 50% starting frequency ~~~~~~~~~~~~~~~~ ##
#!# Not currently used for anything but can be useful for explanatory plots/talks.
chromSNPFreqExpPlot = function(snpdata,chromdata,chrom,alt="Alt",comparison=settings$basename,ymin=0,ymax=1.0,xmin=0,xmax=-1,colbychrom=TRUE,region=FALSE,inverse=FALSE){
  # Plot of the Hits betweeen Chromosomes and Contigs(+) or Chromosomes(-):
  scaling = settings$scaling
  units = paste("(",settings$units,")",sep="")
  if(alt == "Alt"){
    (plotdata = snpdata[snpdata$AltLocus==chrom,])
    plotdata = plotdata[(order(plotdata$AltPos)),]
  }else{
    (plotdata = snpdata[snpdata$Locus==chrom,])
    plotdata = plotdata[(order(plotdata$Pos)),]
  }
  if(inverse){
    plotdata$MajFreq = plotdata$MajFreq - plotdata$MajDiff
  }
  plotdata$MajDiff = plotdata$MajFreq - 0.5
  plotdata$MajProb = 1.0
  if(xmax < xmin){
    xmax = max(chromdata[chromdata$Chrom==chrom,]$End)/scaling
  }
  # Setup Plot
  if(region == FALSE){ region = chrom }
  plot(c(xmin,xmax),c(ymin,ymax),col="red",type="n",ylab="SNP Freq",main=paste(region,comparison),xlab=paste("Chromosome Position",units),xaxs="i")
  plotGridlines(xmin,xmax,ymin,ymax,10000/scaling,50000/scaling,0.1,0.5)
  if(dim(plotdata)[1] < 1){ return(0) }
  # Plot Data
  for(xi in 1:length(plotdata$Pos)){
    pdata = plotdata[xi,]
    # Set plot symbol
    if(pdata$MajDiff > 0){
      ppch=24
    }else if(pdata$MajDiff == 0){
      ppch=23
    }else{
      ppch=25
    }
    # Set plot colours
    if(pdata$MajProb <= settings$pcut){
      p = pdata$MajProb
      i = max(0,(10 + log10(p))/10.0)
      pbg = rgb(1,1,i)
      plwd = 2 - i
    }else{
      pbg = NA
      plwd = 1
    }
    if(colbychrom){
      if(alt == "Alt"){
        pchrom = as.character(pdata$Locus)
        pcol = settings$col[[pchrom]]
      }else{
        pchrom = as.character(pdata$AltLocus)
        pcol = settings$col[[pchrom]]
      }
    }else{
      if(pdata$MajDiff > 0){
        i = min(1,max(0.2,2*pdata$MajDiff))
        pcol = rgb(1,1-i,1-i)
      }else{
        i = max(-1.0,min(-0.2,2*pdata$MajDiff))
        pcol = rgb(1+i,1+i,1)
      }
    }
    # Plot lines and symbols
    if(alt == "Alt"){
      lines(c(pdata$AltPos/scaling,pdata$AltPos/scaling),c(pdata$MajFreq,pdata$MajFreq-pdata$MajDiff),col=pcol,type="l",lwd=plwd)
      points(c(pdata$AltPos/scaling),c(pdata$MajFreq),col=pcol,pch=ppch,bg=pbg)
    }else{
      lines(c(pdata$Pos/scaling,pdata$Pos/scaling),c(pdata$MajFreq,pdata$MajFreq-pdata$MajDiff),col=pcol,type="l",lwd=plwd)
      points(c(pdata$Pos/scaling),c(pdata$MajFreq),col=pcol,pch=ppch,bg=pbg)
    }
  }#lines(plotdata$Pos/scaling,plotdata$ContigNum-plotdata$ChromNum,col="blue",type="l")
}


## ~ Dual panel plot of SNP Frequency and CNV ~~~~~~~~~~~~~~~~ ##
multiSNPCNVPlot = function(pngfile,chrom,alt="Alt",ymin=0,ymax=1.0,xmin=0,xmax=-1,makepng=TRUE){
  if(makepng){
    png(filename=pngfile, width=settings$pngwidth, height=settings$dualpng, units = "px", pointsize=settings$pointsize)
  }
  par(mar=c(2,6,2,1))
  cn = length(settings$freqcomp) + 1
  panels = 1:cn
  layout(matrix(panels,byrow=FALSE,nrow=cn))
  for(ci in 1:(cn-1)){
    if(settings$fdrplot){
      chromSNPFreqPlot(fdrdb[[settings$freqcomp[ci]]],cnvdata,chrom,alt,paste(settings$freqcomp[ci],"FDR"),ymin,ymax,xmin,xmax,colbychrom = TRUE)
    }else{
      chromSNPFreqPlot(combdb[[settings$freqcomp[ci]]],cnvdata,chrom,alt,settings$freqcomp[ci],ymin,ymax,xmin,xmax,colbychrom = TRUE)
    }
  }
  par(mar=c(5,6,1.5,1))
  chromCNVPlot(cnvdata,chrom,xmin=xmin,xmax=xmax) 
  if(makepng){
    dev.off()
  }
}



## ~ Mulitpanel plot of different SNP Frequency changes ~~~~~~~~~~~~~~~~ ##
multiSNPPlot = function(appdata,pngfile,chrom,alt="Alt",ymin=0,ymax=1.0,xmin=0,xmax=-1,makepng=TRUE){
  writeLines("multiSNPPlot")
  plotstatus = "Data plotted."
  # Setup PNG
  if(makepng){
    png(filename=pngfile, width=appdata$pngwidth, height=appdata$pngheight, units = "px", pointsize=appdata$pointsize)
    writeLines(pngfile)
  }
  # Setup Multipanel plot
  par(mar=c(2,6,2,1))
  #i# Original method allowed multiple plots - can bring back: cn = length(settings$freqcomp)
  #!# Add a slider with 1-x comparisons - reveal extra load inputs
  cn = 1  # Temp single comparison
  cx = cn 
  plotft = FALSE
  #if(isFeatureData(alt)){ 
  #  cx = cn + 1 
  #  plotft = TRUE
  #}
  panels = 1:cx
  layout(matrix(panels,byrow=FALSE,nrow=cx))
  for(ci in 1:cn){
    if(ci == cn && plotft == FALSE){
      par(mar=c(5,6,1.5,1))
    }
    if(settings$fdrplot){
      chromSNPFreqPlot(fdrdb[[settings$freqcomp[ci]]],appdata$feattable,chrom,alt,paste(settings$freqcomp[ci],"FDR"),ymin,ymax,xmin,xmax,colbychrom = FALSE)
    }else{
      #chromSNPFreqPlot(combdb[[settings$freqcomp[ci]]],appdata$feattable,chrom,alt,settings$freqcomp[ci],ymin,ymax,xmin,xmax,colbychrom = FALSE)
      plotstatus = chromSNPFreqPlot(appdata,appdata$feattable,chrom,alt,appdata$complabel,ymin,ymax,xmin,xmax,colbychrom = FALSE)
    }
  }
  if(plotft){
    par(mar=c(5,6,1.5,1))
    if(xmax < xmin){
      xmax = max(cnvdata[cnvdata$Chrom==chrom,]$End)/settings$scaling
    }
    featurePlot(alt,chrom,xmin,xmax)
  }  
  if(makepng){
    dev.off()
    return(paste("Saved to ",pngfile,"::",plotstatus))
  }
  return(plotstatus)
}


### ~ Make functions for the different main plot types
pngName = function(appdata){
  pngfile = paste0(appdata$pngpath,"/",appdata$pngbase,".",appdata$chrom,".A-",appdata$parlist,".M-",appdata$strlist,".P-",appdata$poplist,".",appdata$snplist,"-",appdata$efflist)
  if(appdata$freqlist != "Loaded data"){
    pngfile = paste(pngfile,appdata$freqlist,"png",sep=".")
  }
  return(pngfile)
}


#1# Region Zoom plot of a single chromosome
zoomPlot = function(appdata,makepng=FALSE){
  chrom = appdata$chrom
  ext = "png"
  if(settings$fdrplot){ ext = "fdr.png" }
  writeLines(paste(c(chrom,"ZoomPlot")))
  if(settings$chrom %in% appdata$refchr){
    pngfile=pngName(appdata)
    return(multiSNPPlot(appdata,pngfile,chrom,alt="Ref",ymin=settings$ymin,ymax=settings$ymax,xmin=settings$xmin,xmax=settings$xmax,makepng=makepng))
    #pngfile=paste(settings$pngbase,chrom,"cnv",ext,sep=".")
    #multiSNPCNVPlot(pngfile,chrom,alt="Ref",ymin=settings$ymin,ymax=settings$ymax,xmin=settings$xmin,xmax=settings$xmax,makepng=makepng)
  }
  #if(settings$chrom %in% appdata$altchr){
  #  pngfile=paste(settings$pngbase,chrom,ext,sep=".")
  #  multiSNPPlot(appdata,pngfile,chrom,alt="Alt",ymin=settings$ymin,ymax=settings$ymax,xmin=settings$xmin,xmax=settings$xmax,makepng=makepng)
  #  writeLines("multiSNPPlot-Alt")
  #  #pngfile=paste(settings$pngbase,chrom,"cnv",ext,sep=".")
  #  #multiSNPCNVPlot(pngfile,chrom,alt="Alt",ymin=settings$ymin,ymax=settings$ymax,xmin=settings$xmin,xmax=settings$xmax,makepng=makepng)
  #}
}  

#2# Standard Solo plots per chromosome
soloPlots = function(){  
  for(chrom in refchr){
    writeLines(chrom)
    pngfile=paste(settings$pngbase,chrom,"png",sep=".")
    multiSNPPlot(pngfile,chrom,alt="Ref",ymin=settings$ymin,ymax=settings$ymax,xmin=settings$xmin,xmax=settings$xmax)
    #pngfile=paste(settings$pngbase,chrom,"cnv.png",sep=".")
    #multiSNPCNVPlot(pngfile,chrom,alt="Ref",ymin=settings$ymin,ymax=settings$ymax,xmin=settings$xmin,xmax=settings$xmax)
  }
  
  for(chrom in altchr){
    writeLines(chrom)
    pngfile=paste(settings$pngbase,chrom,"png",sep=".")
    multiSNPPlot(pngfile,chrom,alt="Alt",ymin=settings$ymin,ymax=settings$ymax,xmin=settings$xmin,xmax=settings$xmax)
    #pngfile=paste(settings$pngbase,chrom,"cnv.png",sep=".")
    #multiSNPCNVPlot(pngfile,chrom,alt="Alt",ymin=settings$ymin,ymax=settings$ymax,xmin=settings$xmin,xmax=settings$xmax)
    
  }
}

#3# Multiplot of all chromosomes in a single graphic
multiPlotFull = function(makepng=TRUE){
  cn = length(settings$freqcomp)
  pext = "png"
  if(settings$fdrplot){ pext = "fdr.png" }
  for(ci in 1:cn){
    (fcomp = settings$freqcomp[ci])
    comptitle = fcomp
    if(settings$fdrplot){ paste(fcomp,"FDR") }
    if(settings$fdrplot){ plotdata = fdrdb[[fcomp]] }
    else{ plotdata = combdb[[fcomp]] }
    # Reference chromosomes
    alt = "Ref"
    yi = as.integer(length(refchr) / 3)
    if(length(refchr) %% 3 > 0){ yi = yi + 1 }
    height = yi * settings$pngwidth / 6    # Default 400 per plot
    pngfile = paste(settings$pngbase,fcomp,"ref",pext,sep=".")
    if(makepng){
      writeLines(pngfile)
      png(filename=pngfile, width=settings$pngwidth * 1.5, height=height, units = "px", pointsize=settings$pointsize)
    }
    par(mar=c(2,6,2,1))
    par(mar=c(5,6,1.5,1))
    panels = 1:(yi*3)
    layout(matrix(panels,byrow=TRUE,ncol=3))
    for(chrom in refchr[1:(yi*3-4)]){
      chromSNPFreqPlot(plotdata,cnvdata,chrom,alt,comptitle,colbychrom = FALSE)
    }
    par(mar=c(5,6,1,1))
    for(chrom in refchr[(yi*3-3):length(refchr)]){
      chromSNPFreqPlot(plotdata,cnvdata,chrom,alt,comptitle,colbychrom = FALSE)
    }
    if(makepng){
      dev.off()
    }
    
    # Alt chromosomes
    alt = "Alt"
    yi = as.integer(length(altchr) / 3)
    if(length(altchr) %% 3 > 0){ yi = yi + 1 }
    height = yi * settings$pngwidth / 6    # Default 600 per plot
    # SNPFreq plot
    pngfile = paste(settings$pngbase,fcomp,"alt",pext,sep=".")
    if(makepng){
      writeLines(pngfile)
      png(filename=pngfile, width=settings$pngwidth * 1.5, height=height, units = "px", pointsize=settings$pointsize)
    }
    par(mar=c(5,6,1.5,1))
    panels = 1:(yi*3)
    layout(matrix(panels,byrow=TRUE,ncol=3))
    for(chrom in altchr[1:(yi*3-4)]){
      chromSNPFreqPlot(plotdata,cnvdata,chrom,alt,comptitle,colbychrom = FALSE)
    }
    for(chrom in altchr[(yi*3-3):length(altchr)]){
      chromSNPFreqPlot(plotdata,cnvdata,chrom,alt,comptitle,colbychrom = FALSE)
    }
    if(makepng){
      dev.off()
    }
  }  
}


############### ::: TIME SERIES FUNCTIONS ::: ##################
freqTimePlot = function(locus,startPos=1,endPos=0,gloss=TRUE,Title=''){
  if(endPos<1){ 
    endPos = max(snp_table_list_rowsbinded$Pos)
  }
  test = snp_table_list_rowsbinded[snp_table_list_rowsbinded$Locus==locus & snp_table_list_rowsbinded$Pos>=startPos & snp_table_list_rowsbinded$Pos<=endPos,]
  test$SNP = paste(test$Locus,test$Pos,test$AltLocus,test$AltPos,sep=".")
  test$col = "black"
  if(max(test$MajFreq) > 0.95){
    test[test$MajFreq > 0.95,]$col = "red"
  }
  if(min(test$MajFreq) < 0.05){
    test[test$MajFreq < 0.05,]$col = "blue"
  }
  snplist = as.character(levels(as.factor(test$SNP)))
  plot(0,0,type="n",ylim=c(0,1),xlim=c(1,7),ylab="Alternative Allele Freq",xlab="Timepoint",main = Title,xaxt="n")
  axis(1, at=1:7, labels=c("P4","P6","P9","P10","P11","P12","P13"))
  for(snp in snplist){
    onesnp = test[test$SNP==snp,]
    snpcol = onesnp$col[length(onesnp$col)]
    if(length(onesnp$col)==7 | gloss==FALSE){
      lines(onesnp$time_point,onesnp$MajFreq,type="l",col=snpcol)
    }
  }
}




############### ::: REST FUNCTIONS ::: ##################
### Check whether a jobID looks legit and return True or False
isJobID <- function(jobid){
  if(nchar(jobid) != 11){ return(FALSE) }
  if(is.na(as.numeric(jobid))){ return(FALSE) }
  if(is.na(as.integer(substr(jobid,1,6)))){ return(FALSE) }
  if(is.na(as.integer(substr(jobid,7,11)))){ return(FALSE) }
  return(TRUE)
}
### Check whether Job has run
checkJob <- function(jobid,password=""){
  checkurl = paste0(settings$resturl,"check&jobid=",jobid,"&password=",password)
  jobcheck = readLines(checkurl,warn=FALSE)[1]
  if(jobcheck == "Finished"){ return(TRUE) }
  else{ return(jobcheck) }  # Missing/Queued/Running
  #i# Password mismatch: "ERROR: JobID password mismatch! Contact server admin if this is your job and you have forgotten the password used."
} 
### Function for returning the REST keys
getRestKeys <- function(jobid,password=""){
  joburl = paste0(settings$resturl,"retrieve&jobid=",jobid,"&password=",password,"&rest=restkeys")
  return(readLines(joburl,warn=FALSE))
}
### Return an R object with REST output
getRestOutput <- function(jobid,rest,outfmt="",password=""){
  joburl = paste0(settings$resturl,"retrieve&jobid=",jobid,"&password=",password,"&rest=",rest)
  if(outfmt == ""){
    outfmt = "text"
    if(rest %in% settings$csv){ outfmt = "csv"}
    if(rest %in% settings$tsv){ outfmt = "tsv"}
    if(rest == "log"){ outfmt = "log"}
  }
  if(outfmt == "text"){ return(readLines(joburl,warn=FALSE)) }
  if(outfmt == "csv"){ return(read.delim(joburl,header=TRUE,sep=",")) }
  if(outfmt == "tsv"){ return(read.delim(joburl,header=TRUE,sep="\t")) }
  if(outfmt == "log"){ 
    logdata = read.delim(joburl,header=FALSE,sep="\t")
    # Modify Log data
    colnames(logdata) = c("Type","Time","Details")
    logdata$Line = as.integer(rownames(logdata))
    logdata = logdata[,c(4,1,2,3)]
    return(logdata) 
  }
}

############### ::: UPDATE DATA ::: ##################
#!# This is the old function that needs updating with above functions
#i# This function is called by server.R in a reactiveValues() call.
#i# This way, it will be updated whenever the settings that affect it are changed.
setData = function(jobid,prog="retrieve",password="",extra=c(),formats=c()){
  #># jobid = REST server job ID
  #># prog = REST program
  #># password = REST password
  #># extra = Additional 
  restbase = paste("http://rest.slimsuite.unsw.edu.au/",prog,"&jobid=",jobid,"&password=",password,sep="")
  
  # Template SLiMSuite shiny app using example Zen run:
  # http://rest.slimsuite.unsw.edu.au/retrieve&jobid=17092300003&rest=format&password=None&refresh=2
  #i# Standard &rest=X outputs:
  standard = c("status", "version", "ini", "log", "warnings", "errors", "outfmt", "help", "jobid", "intro", "prog")
  #i# Zen-specific &rest=X outputs
  noheadtdt = c("wisdoms")
  headtdt = c()
  noheadcsv = c()
  headcsv = c()
  
  #<# Will return a list of the different REST outputs
  #i# See: https://cran.r-project.org/web/packages/httr/vignettes/quickstart.html
  fulldata = read.delim(paste(restbase,"&rest=full",sep=""),header=FALSE)
  rdata = list(full=fulldata)
  rdata$full = as.character(rdata$full[1:3,1])
  for(rest in c(standard,noheadtdt)){
    #resturl = paste(restbase,"&rest=",rest,sep="")
    rdata[[rest]] <- read.delim(paste(restbase,"&rest=",rest,sep=""),header=FALSE,sep="\t")
  }
  for(rest in c(headtdt)){ rdata[[rest]] <- read.delim(paste(restbase,"&rest=",rest,sep=""),header=TRUE,sep="\t") }
  for(rest in c(noheadcsv)){ rdata[[rest]] <- read.delim(paste(restbase,"&rest=",rest,sep=""),header=FALSE,sep=",") }
  for(rest in c(headcsv)){ rdata[[rest]] <- read.delim(paste(restbase,"&rest=",rest,sep=""),header=TRUE,sep=",") }
  # Modify Log data
  colnames(rdata$log) = c("Type","Time","Details")
  rdata$log$Line = as.integer(rownames(rdata$log))
  rdata$log = rdata$log[,c(4,1,2,3)]
  # Return data list
  return(rdata)
}
#cat(as.character(pep[1:10,1]))

############### ::: SHINY CODE INFO ::: ##################
#i# This section contains some information comments that can be deleted in actual Apps.

# Template SLiMSuite shiny app using example Zen run:
# http://rest.slimsuite.unsw.edu.au/retrieve&jobid=17092300003&rest=format&password=None&refresh=2
#i# Standard &rest=X outputs:
# - status, version, ini, log, warnings
#i# Zen-specific &rest=X outputs
# - wisdoms

# # status: ./17092300003.status
# JobID 17092300003 (zen) Finished.
# IP:58.178.74.181
# No queue.
# zen&wisdoms=10
# Run Started: 2017-09-23 13:37:10; PID=5456
# Run finished: 2017-09-23 13:37:11

#i# Output options
#i# Render text inside <pre> html block:
# verbatimTextOutput(outputId, placeholder = FALSE) 

# htmlOutput (uiOutput) = Create an HTML output element
# plotOutput (imageOutput) = Create an plot or image output element
# outputOptions = Set options for an output object.
# tableOutput (dataTableOutput) = Create a table output element
# textOutput = Create a text output element
# verbatimTextOutput = Create a verbatim text output element
# downloadButton (downloadLink) = Create a download button or link
# Progress = Reporting progress (object-oriented API)
# withProgress (setProgress, incProgress) = Reporting progress (functional API)
# modalDialog = Create a modal dialog UI
# urlModal = Generate a modal dialog that displays a URL
# showModal (removeModal) = Show or remove a modal dialog
# showNotification (removeNotification) = Show or remove a notification

# renderPlot = Plot Output
# renderText = Text Output
# renderPrint = Printable Output
# renderDataTable = Table output with the JavaScript library DataTables
# renderImage = Image file output
# renderTable = Table Output
# renderUI = UI Output
# downloadHandler = File Downloads

#i# UI Inputs:
# actionButton (actionLink) = Action button/link
# checkboxGroupInput = Checkbox Group Input Control
# checkboxInput = Checkbox Input Control
# dateInput = Create date input
# dateRangeInput = Create date range input
# fileInput = File Upload Control
# numericInput = Create a numeric input control
# radioButtons = Create radio buttons
# selectInput (selectizeInput) = Create a select list input control
# sliderInput (animationOptions) = Slider Input Widget
# submitButton = Create a submit button
# textInput = Create a text input control
# textAreaInput = Create a textarea input control
# passwordInput = Create a password input control
# modalButton = Create a button for a modal dialog
# updateActionButton = Change the label or icon of an action button on the client
# updateCheckboxGroupInput = Change the value of a checkbox group input on the client
# updateCheckboxInput = Change the value of a checkbox input on the client
# updateDateInput = Change the value of a date input on the client
# updateDateRangeInput = Change the start and end values of a date range input on the client
# updateNumericInput = Change the value of a number input on the client
# updateRadioButtons = Change the value of a radio input on the client
# updateSelectInput (updateSelectizeInput) = Change the value of a select input on the client
# updateSliderInput = Change the value of a slider input on the client
# updateTabsetPanel (updateNavbarPage, updateNavlistPanel) = Change the selected tab on the client
# insertTab (prependTab, appendTab, removeTab) = Dynamically insert/remove a tabPanel
# showTab (hideTab) = Dynamically hide/show a tabPanel
# updateTextInput = Change the value of a text input on the client
# updateTextAreaInput = Change the value of a textarea input on the client
# updateQueryString = Update URL in browser's location bar
# getQueryString (getUrlHash) = Get the query string / hash component from the URL
