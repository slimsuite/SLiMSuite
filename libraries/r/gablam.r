#########################################
### GABLAM.R - SCRIPT FOR GABLAM V2.9 ###
### AUTHOR: Dr Richard Edwards 2011	~ ###
### CONTACT: r.edwards@soton.ac.uk ~~ ###
### Last Edit: 25/7/12 ~~~~~~~~~~~~~~ ###
#########################################

############# ::: GENERAL SETUP ::: ##################

useCairo = FALSE	               # Not for bioinf.soton.ac.uk
if(useCairo){
	library(Cairo)    
}

## Setup Path to general R Files
rdir = "/home/re1u06/researchfiles/SBSBINF/Tools/svn/libraries/r/"
#rdir = "/home/bioinf/bioware/Tools/libraries/r/"
#rdir = "T:\\Tools\\svn\\libraries\\r\\"
rdir = "/Users/redwards/Dropbox/_Repository_/slimsuite/libraries/r/"


rjesource = function(rfile){
    source(paste(rdir,rfile,sep=""))
}
rjesource("rje_col.r")
rjesource("rje_misc.r")

### ~ COMMANDLINE ARGUMENTS ~ ###
args <- commandArgs(TRUE)
(basefile = args[1])
if(length(args) > 1){
	(minfraglen = as.integer(args[2]))
}else{
	minfraglen = 0
}


############# ::: LOAD DATA FOR ANALYSIS ::: #########
gabfile = paste(basefile,"gablam.tdt",sep=".")
gabdata = read.table(gabfile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
summary(gabdata)
locfile = paste(basefile,"local.tdt",sep=".")
locdata = read.table(locfile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
summary(locdata)

############# ::: PLOTTING FUNCTION ::: #########
dotplot = function(qry,hit,makepng=TRUE,pcol=rep(soton$col[1],100),minfrag=0){
	### ~ Setup Plot ~ ###
	pngwidth = 1600
	outfile = paste(basefile,".DotPlots//",qry,".",hit,".dot.png",sep="")
	if (makepng){
		if (useCairo){
			CairoPNG(filename=outfile, width = pngwidth, height = pngwidth, pointsize=24)
		}else{
			png(filename = outfile, width = pngwidth, height = pngwidth, units = "px", pointsize=24)
		}
	}
	(gdata = gabdata[gabdata$Qry==qry & gabdata$Hit==hit,])
	(qlen = gdata$QryLen)
	(hlen = gdata$HitLen)
	plot(c(1,qlen),c(1,hlen),type="n",xlab=qry,ylab=hit,axes=TRUE,ann=TRUE,mar=c(0,1,4,1))
	title(main=paste(hit," vs ",qry," (Minfrag=",minfrag,"bp)",sep=""))
	### ~ Plot local hit data ~ ###
	ldata = locdata[locdata$Qry==qry & locdata$Hit==hit & locdata$Length >= minfrag,]
	ldata$Perc = ldata$Identity / ldata$Length
	for(i in 1:length(ldata$Perc)){
		(tran = as.integer(100*ldata$Perc[i]) - 1)
		#(ldata[i,])
		#icol = (paste(pcol,tran,sep=""))
		icol = pcol[tran+1]
		lines(c(ldata$QryStart[i],ldata$QryEnd[i]),c(ldata$SbjStart[i],ldata$SbjEnd[i]),type="l",lwd=2,col=icol)
		lines(c(ldata$QryStart[i],ldata$QryEnd[i]),c(ldata$SbjStart[i],ldata$SbjEnd[i]),type="p",lwd=2,col=icol)
	}
	if (makepng){
		dev.off()
	}
}


############# ::: PERFORM PLOTS ::: #########
pcol = c(rep("red",89),rep("blue",10),"black")
(plotnum = length(gabdata$Qry))
for(g in 1:plotnum){
	qry = gabdata[g,]$Qry
	hit = gabdata[g,]$Hit
	dotplot(qry,hit,TRUE,pcol,minfraglen)
}

#dotplot(qry,hit,FALSE,pcol,1000)