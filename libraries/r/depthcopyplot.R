########################################################
### DepthKopy plot functions                   ~~~~~ ###
### VERSION: 0.1.0                             ~~~~~ ###
### LAST EDIT: 01/11/22                        ~~~~~ ###
### AUTHORS: Richard Edwards 2022              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for DepthKopy copy plotting from compiled data

####################################### ::: HISTORY ::: ############################################
# v0.0.0 : Initial dev version based on depthcopy.R
# v0.1.0 : Initial working version with plotting from DepthKopy Excel output.
version = "v0.1.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript depthcopyplot.R basefile=FILE [scdepth=NUM] [xlsx=FILE] [xsheets=LIST] [reghead=LIST] [pngdir=PATH] [cnmax=INT] [sigdif=T/F] [rdir=PATH] 
# : basefile=PREFIX = Basefile prefix. Will look for $BASEFILE.xlsx to load data if not given.
# : scdepth=NUM = Single copy read depth to use in place of BUSCO-derived depth
# : xlsx=FILE = Full XLSX results table to load. (Use if different to $BASEFILE.xlsx)
# : xsheets=LIST = Optional list of xlsx sheets to include in plot
# : reghead=LIST = Optional SeqName,Start,End or GFF feature type. Defaults to SeqName,Start,End or "gene"
# : pngdir=PATH = output path for PNG graphics
# : cnmax=INT = maximum CN value for output graphics
# : sigdif=T/F = Whether to calculate statistical differences between data [False]
# : rdir=PATH = Path to other R scripts. (Will work out from Rscript call if not provided)

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# os.popen('Rscript {0}depthcopy.R {1}'.format(rdir, optionstr)).readlines()

####################################### ::: OUTPUTS ::: ############################################
#!# List the outputs here

####################################### ::: TODO ::: ############################################
#!# Tidy up input

####################################### ::: SETUP ::: ############################################
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(scdepth=0,debug=FALSE,cnmax=4,xlsx="",xsheets=vector(),
                seqstats=FALSE,sigdif=FALSE,
                pngwidth=1200,pngheight=900,pointsize=16,pngdir="",
                basefile="depthcopy",rdir="",
                reghead="SeqName,Start,End",outlog=stdout())
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
for(cmd in c("pngwidth","pngheight","pointsize","cnmax")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
for(cmd in c("scdepth")){
  settings[[cmd]] = as.numeric(settings[[cmd]])
}
for(cmd in c("reghead","xsheets")){
  if(length(settings[[cmd]]) > 0){
    settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
  }else{
    settings[[cmd]] <- defaults[[cmd]]
  }
}
for(cmd in c("debug","sigdif")){
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
logWrite(paste("#RCODE DepthCopyPlot.R:",version))
for(cmd in names(settings)[order(names(settings))]){
  logWrite(paste("CMD:",cmd,"=",paste0(settings[[cmd]],collapse=",")))
}
#x#settings$buscofile = settings$busco
if(settings$pngdir == ""){
  settings$pngdir <- paste0(settings$basefile,".plots")
}
dir.create(settings$pngdir, showWarnings = FALSE)

### ~ Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(! "tidyverse" %in% installed.packages()[,"Package"]){
  install.packages("tidyverse")
}
library(tidyverse)
library(RColorBrewer)
settings$ggstatsplot = "ggstatsplot" %in% installed.packages()[,"Package"]
if(settings$ggstatsplot){
  library(ggstatsplot)
}
settings$readxl = "readxl" %in% installed.packages()[,"Package"]
if(settings$readxl){
  library(readxl)
}else{
  logWrite("#XLXS Install readxl package for compiled Excel file input")
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
#sfile <- paste0(settings$rdir,"/rje_load.R")
#logWrite(sfile)
#source(sfile)

####################################### ::: FUNCTIONS ::: ############################################

##### ======================== Loading data functions ======================== #####

### ~ Load Excel File into full data frame ~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Loads Excel file and generates the plotting database
loadPlotXL <- function(filename){
  # Setup plotdb
  newplotdb <- tibble(Dataset=c(),SeqName=c(),Start=c(),End=c(),MeanX=c(),MedX=c(),ModeX=c(),DensX=c(),CN=c(),DensK=c(),SelfK=c(),HomPC=c())
  reghead <- settings$reghead
  # Load and Compile Excel workbook
  for(xsheet in excel_sheets(filename)){
    #print(xsheet)
    xldb <- read_excel(filename, sheet=xsheet)

    if(colnames(xldb)[1] == "BuscoID"){
      xldb <- xldb %>% rename(SeqName=Contig)
    }
    if(sum(settings$reghead %in% colnames(xldb)) == 3){
      xldb <- xldb %>% rename(SeqName=reghead[1],Start=reghead[2],End=reghead[3])
    }
    #i# Assume first three fields are SeqName, Start, End
    if(sum(c("SeqName","Start","End") %in% colnames(xldb)) != 3){
      xldb <- xldb %>% rename(SeqName=colnames(xldb)[1],Start=colnames(xldb)[2],End=colnames(xldb)[3])
    }
      
    #print(colnames(xldb))
    
    newplotdb <- bind_rows(newplotdb, 
                           xldb %>% 
                             mutate(Dataset=xsheet) %>%
                             select(Dataset,SeqName,Start,End,MeanX,MedX,ModeX,DensX,CN,DensK,SelfK,HomPC))
    
    logWrite(paste("#XLXS Loaded",nrow(xldb),"from sheet",xsheet,"->",nrow(newplotdb),"full data points."))
    
  }
  return(newplotdb)
  
}

### ~ Tidy tables for output ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Sets dp etc. for output fields in a table then returns for output
tidyTable <- function(regdb){
  ifields = c("Start","End",reghead[2],reghead[3],"ModeX")
  dp2fields = c("MeanX","MedX","DensX","DensK","HomPC")
  dp3fields = c("SelfK","CN","CIsyst","CIrand")
  #i# Integers
  for(colname in ifields){
    if(colname %in% colnames(regdb)){
      regdb[[colname]] = as.integer(regdb[[colname]])
    }
  }
  #i# 2 d.p.
  for(colname in dp2fields){
    if(colname %in% colnames(regdb)){
      regdb[[colname]] = round(regdb[[colname]],2)
    }
  }
  #i# 3 d.p.
  for(colname in dp3fields){
    if(colname %in% colnames(regdb)){
      regdb[[colname]] = round(regdb[[colname]],2)
    }
  }
  return(regdb)
}

##### ======================== Plotting functions ======================== #####

### ~ Generate Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Violin plot function
vioPlot <- function(plotdb,plotfield){
  #i# Violin plots by Method
  #!# Add labels of n=X for each Dataset
  labels = levels(factor(plotdb$Dataset))
  labels = paste("n =",table(factor(plotdb$Dataset)))
  nx = 1:length(labels)
  ny = rep(0,length(labels))
  p <- ggplot(plotdb, aes(factor(Dataset), .data[[plotfield]]), palette = palette) + 
    aes(fill = factor(Dataset)) +
    geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) + 
    xlab("Dataset") + ylab(plotfield) +
    labs(fill="Dataset") +
    #geom_text(x=nx, y=ny, label=labels, vjust=1) +
    theme_bw(base_size = settings$pointsize)
  
  if(settings$ggstatsplot){
    splotdb = plotdb
    splotdb$plotfield = splotdb[[plotfield]]
    p <- ggbetweenstats(
      data = splotdb,
      results.subtitle = length(labels) <= 4 & settings$sigdif,
      pairwise.comparisons = settings$sigdif,
      x = Dataset,
      y = plotfield
    ) + ylim(0,median(splotdb$plotfield)*settings$cnmax) +
      #geom_hline(yintercept=1, color="black") +
      labs(
        #title = "DepthSizer accuracy versus reference assembly size",
        x = "Dataset",
        y = plotfield,
        caption = NULL
      ) +
      theme(text = element_text(size=settings$pointsize))
    
  }
  
  nlevel = length(levels(factor(plotdb$Dataset)))
  if(nlevel > 8){
    p <- p + scale_colour_manual(values=rep(brewer.pal(12,"Paired"),times=as.integer(nlevel/12)+1)[1:nlevel])
  }
  
  if("Other" %in% levels(factor(plotdb$Dataset))){
    p <- p + theme(axis.text.x = element_text(angle = 90))
  }
  
  if(plotfield %in% c("CN","SelfK")){
    p <- p + geom_hline(yintercept=1, color="steelblue", linetype="dashed")
  }
  if(plotfield %in% c("MeanX","MedX","ModeX","DensX")){
    p <- p + geom_hline(yintercept=scdepth, color="steelblue", linetype="dashed")
  }
  
  return(p)
}

### ~ Depth / Kmer plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
hexPlot <- function(plotdb,plotfield){
  p = ggplot(plotdb, aes(x=CN, y=.data[[plotfield]]) ) +
    geom_hex(bins = 100) + #scale_x_log10() + scale_y_log10() +
    xlab("CN") + ylab(plotfield) +
    #coord_trans(x = "log2", y = "log2") +
    scale_x_continuous(trans = "log2") + 
    scale_fill_continuous(type = "viridis") +
    geom_vline(xintercept=1, color="black") +
    # geom_hline(yintercept=0.5, linetype=2, color="#d29f23") +
    # geom_vline(xintercept=0.5, linetype=2, color="#d29f23") +
    # geom_hline(yintercept=2, linetype=2, color="#d29f23") +
    # geom_vline(xintercept=2, linetype=2, color="#d29f23") +
    theme_bw(base_size = settings$pointsize)
  
  if(! plotfield %in% c("HomPC")){
    p = p +  scale_y_continuous(trans = "log2") 
  }
  
  return(p)
}

################################## ::: RUN CODE ::: #######################################

##### ======================== Report key inputs ======================== #####

### ~ Key inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# SC Depth
scdepth = as.numeric(settings$scdepth)
if(scdepth > 0){
  logWrite(paste("SC Depth:",round(scdepth,2)))
}else{
  logWrite("No scdepth=NUM BUSCO CN depth given. Will not be marked on X depth plots.")
  scdepth <- 0
  #?# Why not just calculate from CN field?
}
#i# Excel file of all the depth and CN calculations
xfile <- settings$xlsx
if(xfile == ""){
  xfile <- paste0(settings$basefile,".xlsx")
}
if(file.exists(xfile)){
  logWrite(paste("Excel Full File:",xfile))
}else{
  logWrite(paste("Excel Full File missing:",xfile))
  stop("No Excel file to load!")
}
#i# Region file headers
reghead <- settings$reghead
logWrite(paste("Region Headers:",paste0(reghead,collapse=", ")))
#i# End of checks/setups
logWrite('#RCODE Setup complete.')


##### ======================== Load data ======================== #####

### ~ Generate new Plot database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
newplotdb <- loadPlotXL(xfile)
#i# Restrict to the sheets to analysis
if(length(settings$xsheets) > 0){
  xlevels <- settings$xsheets
}else{
  xlevels <- excel_sheets(xfile)
}
xlevels <- xlevels[xlevels %in% newplotdb$Dataset]
if(length(xlevels) > 0){
  logWrite(paste('#DSETS',length(xlevels),"dataset factors to plot:",paste(xlevels,collapse=", ")))
}else{
  stop("No dataset types retained to plot.")
}
#i# Filter and convert to ordered Dataset factors
newplotdb <- newplotdb[newplotdb$Dataset %in% xlevels,]
newplotdb$Dataset <- ordered(newplotdb$Dataset,levels=xlevels)

#!# Plotting code from depthcopy.R
plotdb <- newplotdb
outfile = paste0(settings$basefile,".fulltable.tsv")
logWrite(paste("#SAVE Compiled CN statistics for plotting output to",outfile))
write.table(tidyTable(plotdb),outfile,sep="\t",quote=FALSE,row.names=FALSE)


##### ======================== Plotting ======================== #####

#i# Generate PNG of violin plots
plotfields = c("MeanX","MedX","ModeX","DensX","CN")
for(kfield in plotfields){
  nax <- sum(is.na(plotdb[[kfield]]))
  if(nax > 0){
    plotdb[[kfield]][is.na(plotdb[[kfield]])] <- 0
    logWrite(paste("#WARN",nax,kfield,"NA entries replaced with 0."))
  }
}
for(kfield in c("DensK","SelfK","HomPC")){
  nax <- sum(is.na(plotdb[[kfield]]))
  if(nax > 0){
    plotdb[[kfield]][is.na(plotdb[[kfield]])] <- 0
    logWrite(paste("#WARN",nax,kfield,"NA entries replaced with 0."))
  }
  if(max(plotdb[[kfield]]) > 0){
    plotfields <- c(plotfields, kfield)
  }
}

if(nrow(plotdb) > 0){
  for(plotfield in plotfields){
    if(! plotfield %in% colnames(plotdb)){
      next
    }
    if(sum(is.na(plotdb[[plotfield]]))>=1){
      logWrite(paste("Missing",plotfield,"data: no violin plot."))
      next
    }
    if(max(plotdb[[plotfield]]) > 0){
      plt <- vioPlot(plotdb,plotfield)
      pdffile = paste(settings$basefile,plotfield,"pdf",sep=".")
      pdffile = file.path(settings$pngdir,pdffile)
      ggsave(
        filename = pdffile,
        plot = plt + theme(text = element_text(size=settings$pointsize)),
        device = "pdf",
        width = settings$pngwidth/30,
        height = settings$pngheight/30,
        units = "cm"
      )
      logWrite(paste(plotfield,"violin plot output to:",pdffile))
      pngfile = paste(settings$basefile,plotfield,"png",sep=".")
      pngfile = file.path(settings$pngdir,pngfile)
      png(pngfile,width=settings$pngwidth,height=settings$pngheight,pointsize=settings$pointsize*2)
      print(plt)
      dev.off()
      logWrite(paste(plotfield,"violin plot output to:",pngfile))
    }else{
      logWrite(paste("No",plotfield,"data for violin plot."))
    }
  }
}

#i# Generate PNG of hex plots
hexdir = file.path(settings$pngdir,'hexplots')
dir.create(hexdir, showWarnings = FALSE)
for(plotfield in plotfields[plotfields != "CN"]){
  for(dset in unique(plotdb$Dataset)){
    datdb = plotdb %>% filter(Dataset==dset)
    if(nrow(datdb) > 0){
      plt <- hexPlot(datdb,plotfield)
      pngfile = paste(settings$basefile,dset,plotfield,"png",sep=".")
      pngfile = file.path(hexdir,pngfile)
      png(pngfile,width=settings$pngwidth,height=settings$pngheight,pointsize=settings$pointsize*2)
      print(plt)
      dev.off()
      logWrite(paste(dset,plotfield,"hex plot output to:",pngfile))
    }
  }
}
logWrite(paste0('#PLOTS PNG graphics output to ',settings$pngdir,"/",settings$basefile,".*"))
if(settings$ggstatsplot){
  logWrite('#CITE If using violin plots, please cite:  Patil, I. (2021). J. Open Source Software, 6(61), 3167, doi:10.21105/joss.03167')
}

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE DepthCopyPlot.R finished.")


