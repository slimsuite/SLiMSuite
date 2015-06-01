### BUDAPEST TREE R SCRIPT ###

#rdir = "/home/re1u06/researchfiles/SBSBINF/Tools/_Serpentry/"
#rdir = "/home/redwards/Serpentry/"
#rdir = "/data/ben/Serpentry/"
rdir = "/home/re1u06/researchfiles/SBSBINF/Tools/svn/libraries/r/"
rdir = "/Users/redwards/Dropbox/_Repository_/slimsuite/libraries/r/"

rjesource = function(rfile){
    source(paste(rdir,rfile,sep=""))
}
rjesource("rje_col.r")
rjesource("rje_misc.r")

### ~ COMMANDLINE ARGUMENTS ~ ###
args <- commandArgs(TRUE)
(file = args[1])
(outfile = args[2])		# "KCMA1.png"
if (length(args) > 2) { treetitle = args[3]; }


### ~ SETUP ~ ###
#library(Cairo)                   # Not for bioinf.soton.ac.uk
tree = read.csv(file)
ynum = length(tree$xpos)
x = 1 / (2*max(tree$xpos))
y = 1 / (ynum+1)

### Families and colours ###
fcolist = c(1:4,8:21)
tree$fcol = "black"       # soton$col[1]	#"black"
flvl = levels(as.factor(tree[tree$family != "bud",]$family))
#flvl = levels(as.factor(tree$family))
if(length(flvl) > 0){
                for (f in 1:length(flvl)){ tree[tree$family == flvl[f],]$fcol = soton$col[fcolist[f]] }
}
#~special~#
#if (length(tree[tree$family == "EHUX",]) > 0){
#   tree[tree$family == "EHUX",]$fcol = soton$col[5]
#}

for(subfam in levels(as.factor(tree$family))){
    for(qryfam in c("bud","EHUX","qry","fdr","cand")){
        if(subfam == qryfam){ tree[tree$family == qryfam,]$fcol = soton$col[5] }
    }
}

### ~ Setup Plot ~ ###
pngwidth = 1600
yex = 100 #38
ypx = 20 #50
mex = 2  #1
if (length(args) > 3) { pngwidth = as.integer(args[4]); }
png(filename = outfile, width = pngwidth, height = max(300,(ynum+2)*ypx), units = "px", pointsize=12)
#CairoPNG(filename=outfile, width = pngwidth, height = max(300,(ynum+2)*ypx), pointsize=12)
plot(0:1,0:1,type="n",axes=FALSE,ann=FALSE,mar=c(0,1,4,1))
if (length(args) > 2) { title(main=treetitle); }

### ~ Draw Tree ~ ###
for (i in 1:ynum){
	data = tree[i,]
	lines(c(data$xpos*x,data$ancx*x),c(1-data$ypos*y,1-data$ypos*y),col=data$fcol)
	lines(c(data$ancx*x,data$ancx*x),c(1-data$ypos*y,1-data$ancy*y),col=data$fcol)
}
### ~ Add Text ~ ###
for (i in 1:ynum){
	data = tree[i,]
	if (data$nodenum <= ((ynum+1) / 2)){
		#text((data$xpos)*x+0.01,1-data$ypos*y,data$name,adj=c(0,0.5),cex=min(mex,(yex+2)/ynum),col=data$fcol)
		text((data$xpos)*x+0.01,1-data$ypos*y,data$name,adj=c(0,0.5),col=data$fcol)
	}else{
		#text((data$xpos)*x-0.01,1+(0.45/ynum)-data$ypos*y,data$boot,adj=c(1,0),cex=min(mex,yex/ynum),col=soton$col[5])
		text((data$xpos)*x-0.01,1+(0.45/ynum)-data$ypos*y,data$boot,adj=c(1,0),cex=0.8,col=soton$col[5])
	}
}
### ~ Add scale ~ ###
lines(c(0,0.1*x),c(0,0),col=soton$col[7])
lines(c(0,0),c(0,-0.005),col=soton$col[7])
lines(c(0.1*x,0.1*x),c(0,-0.005),col=soton$col[7])
#text(0,-0.01,"0",adj=c(0.5,1),cex=min(mex,yex/ynum),col=soton$col[7],xpd=NA)
#text(0.1*x,-0.01,"0.1",adj=c(0.5,1),cex=min(mex,yex/ynum),col=soton$col[7],xpd=NA)
text(0,-0.01,"0",adj=c(0.5,1),col=soton$col[7],xpd=NA)
text(0.1*x,-0.01,"0.1",adj=c(0.5,1),col=soton$col[7],xpd=NA)

dev.off()
