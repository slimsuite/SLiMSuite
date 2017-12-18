######################################################
### GENERAL USAGE FUNCTIONS FOR RJE R SCRIPTS ~~~~ ###
### AUTHOR: Dr Richard Edwards 2008 ~~~~~~~~~~~~~~ ###
### CONTACT: r.edwards@soton.ac.uk ~~~~~~~~~~~~~~~ ###
######################################################

############# ::: DATA FUNCTIONS ::: ################
loadData = function(datadb=list(),baselist=c(),extlist=c(),inpath=""){
  for(fbase in baselist){
    datadb[[fbase]] = list()
    for(fext in extlist){
      (tdtfile = paste(inpath,fbase,".",fext,".tdt",sep=""))
      if(file.exists(tdtfile)){
        datadb[[fbase]][[fext]] = read.table( tdtfile, header = T, stringsAsFactors = F, sep = '\t', quote = '', comment.char="")
        summary(datadb[[fbase]][[fext]])
      }
    }
  }
  return(datadb)  
}

getData = function(datadb,fbase,fext){
  if(length(datadb[[fbase]][[fext]])>0){
    return(datadb[[fbase]][[fext]])
  }else{
    return(data.frame())
  }
}

checkData = function(baselist=c(),extlist=c(),inpath=""){
  checkdb = list()
  for(fbase in baselist){
    checkdb[[fbase]] = list()
    for(fext in extlist){
      (tdtfile = paste(inpath,fbase,".",fext,".tdt",sep=""))
      checkdb[[fbase]][[fext]] = c(tdtfile,file.exists(tdtfile))
    }
  }
  return(checkdb)  
}

############# ::: GENERAL FUNCTIONS ::: ################

### ~ Adding leading zeros to numbers and returning as string ~ ###
preZero = function(n,digits=5){
	nstr = as.character(n)
	while (nchar(nstr)<digits){ nstr = paste("0",nstr,sep="") }
	return(nstr)
}



############# ::: PLOTTING FUNCTIONS ::: ################

### ~ Plot Gridlines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
plotGridlines = function(xmin=0,xmax=0,ymin=0,ymax=0,xminor=0,xmajor=0,yminor=0,ymajor=0,xzero=TRUE,yzero=TRUE,xcap=FALSE,ycap=FALSE){
  # Note that xmin and ymin should only be given if negative
  if(xmin < 0){
    if(xminor > 0){ abline(v=c(1:(xmin/xminor))*xminor,lty="dotted",col="grey") }
    if(xmajor > 0){ abline(v=c(1:(xmin/xmajor))*xmajor,col="grey") }
  }
  if(xmax > 0){
    if(xminor > 0){ abline(v=c(1:(xmax/xminor))*xminor,lty="dotted",col="grey") }
    if(xmajor > 0){ abline(v=c(1:(xmax/xmajor))*xmajor,col="grey") }
  }
  if(xzero == TRUE){ abline(v=0,col="grey") }
  if(xcap == TRUE){ abline(v=xmax,col="grey") }
  if(ymin < 0){
    if(yminor > 0){ abline(h=c(1:(ymin/yminor))*yminor,lty="dotted",col="grey") }
    if(ymajor > 0){ abline(h=c(1:(ymin/ymajor))*ymajor,col="grey") }
  }
  if(ymax > 0){
    if(yminor > 0){ abline(h=c(1:(ymax/yminor))*yminor,lty="dotted",col="grey") }
    if(ymajor > 0){ abline(h=c(1:(ymax/ymajor))*ymajor,col="grey") }
  }
  if(yzero == TRUE){ abline(h=0,col="grey") }
  if(ycap == TRUE){ abline(h=ymax,col="grey") }
}
