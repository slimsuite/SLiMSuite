######################################################
### GENERIC COLOUR PALETTE SETUP SCRIPT FOR RJE.R  ###
### AUTHOR: Dr Richard Edwards 2008 ~~~~~~~~~~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au ~~~~~~~~~~~~~~~ ###
######################################################

## ~ University of Southampton Colour Palette ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
soton = list(col = c("#014359",	# Marine Blue
	"#007275",	# Turquoise
	"#0A96A9",	# Light blue
	"#323D43",	# Dark grey
	"#979E45",	# Gold
	"#BBBBBB","#9BA3A6",	# Metallic ~ 1y End ~
	"#653A28",	# Dark Red
	"#531F44",	# Maroon
	"#A67891",	# Dusky Pink
	"#B2699F",	# Pink
	"#CCDAEA",	# Pale light blue ~ 2y1 End ~
	"#8A412B",	# Brown
	"#AB1210",	# Brick red
	"#F00F2C",	# Red
	"#FE3E14",	# Orange
	"#FFB300",	# Yellow ~ 2y2 End ~
	"#4F5A20",	# Green
	"#91BA91",	# Pale green
	"#BDB68A",	# Pale brown
	"#8F9E94"))	# Green/grey ~ End
soton$primary = soton$col[1:7]
soton$secondary = soton$col[8:length(soton$col)]
soton$sec1 = soton$col[8:12]
soton$sec2 = soton$col[13:17]
soton$sec3 = soton$col[18:21]

## ~ Amino acid colour palette ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
acol = list()
for (aa in c("I","L","V","IFLMV","IFLM","IFMV","IFLV","FLMV","ILMV","IL","IV","ILV","ILF")){ acol[aa] = soton$col[19] }
for (aa in c("A","M")){ acol[aa] = soton$col[19] }
for (aa in c("AG","GS","AS","AGS")){ acol[aa] = soton$col[21] }
for (aa in c("K","R","KR")){ acol[aa] = soton$col[3] }
for (aa in c("D","E","DE")){ acol[aa] = soton$col[15] }
for (aa in c("S","T","ST")){ acol[aa] = soton$col[14] }
for (aa in c("C")){ acol[aa] = soton$col[9] }
for (aa in c("P")){ acol[aa] = soton$col[17] }
for (aa in c("F","Y","W","FY","FW","WY","FWY")){ acol[aa] = soton$col[18] }
for (aa in c("G")){ acol[aa] = soton$col[5] }
for (aa in c("H","HK","HR","HKR")){ acol[aa] = soton$col[2] }
for (aa in c("HY","FH","FHY")){ acol[aa] = soton$col[1] }
for (aa in c("Q","N")){ acol[aa] = soton$col[1] }
for (aa in c("X","x",".")){ acol[aa] = soton$col[6] }
for (aa in c("0","1","2","3","4","5","6","7","8","9","[","]",",","{","}")){ acol[aa] = "white" }
for (aa in c("-")){ acol[aa] = "black" }

## ~ PPI Network Colours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
ppcol = c("black",soton$col[1:11],soton$col[13:21])

## ~ Colour Test Function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
colTest = function(){
	plot(0:1,0:1,type="n",axes=FALSE,ann=FALSE,mar=c(0.1,0.5,0.5,0.1))
	x = 1.0 / length(soton$col)
	for (i in 1:length(soton$col)){
		rect((i-1)*x,0.5,i*x,1,col=soton$col[i])
		text((i-0.5)*x,0.75,i,adj=c(0.5,0.5),font=2)#,col=soton$col[17],cex=2)
	}
	x = 1.0 / 21
	i = 1
	for (aa in c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y")){
		rect((i-1)*x,0,i*x,0.5,col=acol[[aa]])
		text((i-0.5)*x,0.25,aa,adj=c(0.5,0.5))
		i = i + 1
	}
}
#colTest()

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


