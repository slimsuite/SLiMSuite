######################################################
### COMPOSITE IMAGE SCRIPT TO ACCOMPANY SLIMJIM.PY ###
### AUTHOR: Dr Richard Edwards 2008	~~~~~~~~~~~~~~ ###
### CONTACT: r.edwards@soton.ac.uk ~~~~~~~~~~~~~~~ ###
######################################################

############# ::: GENERAL SETUP ::: #####################################################################################

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

#ppcol = c("black",soton$col[1:11],soton$col[13:16],soton$col[18:21])
ppcol = c("black",soton$col[1:11],soton$col[13:21])


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

############# ::: CREATE FUNCTIONS ::: ##################################################################################

preZero = function(n,digits=5){
	nstr = as.character(n)
	while (nchar(nstr)<digits){ nstr = paste("0",nstr,sep="") }
	return(nstr)
}

## ~ Generate sequence logo from profile data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
seqLogo = function(infile,fex=1){
	## ~ Read Data and Setup ~ ##
	motif = read.table(infile,sep="\t",header=TRUE)
	aas = colnames(motif)[-1]
	x = 1.0 / length(motif[,1])
	## ~ Create empty plot ~ ##
	old_mar = par()$mar
	par(mar=c(0,2,2,0))
	#!# png(filename = outfile, width=50*length(motif[,1]), height=500, units = "px")	#, pointsize=25)
	plot(0:1,0:1,type="n",axes=FALSE,ann=FALSE,mar=c(0.1,0.5,0.5,0.1))
	axis(side=2,lwd=2,cex=fex)
	## ~ Plot Profile ~ ##
	ix = 0.0
	for (i in 1:length(motif[,1])){
		## Normalise data ##
		isum = sum(motif[i,-1])
		for (a in aas){ motif[i,a] = motif[i,a] / isum }
		## Plot ##
		iy = 0.0
		for (a in order(rank(motif[i,]))){
			aa = colnames(motif)[a]
			if (aa == "Pos"){
				text(ix+(x/2),1.05,motif[i,1],adj=c(0.5,0.5),xpd=TRUE,font=2,family="mono",cex=fex*min(1,4/nchar(motif[i,1])))
			}else{
				f = motif[i,aa]
				rect(ix,iy,ix+x,iy+f,col=acol[[aa]],lwd=2)
				if (f > 0.05){ text(ix+(x/2),iy+(f/2),aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex) }
				iy = iy + f
			}
		}
		ix = ix + x
	}
	#!# dev.off()
	par(mar=old_mar)
	length(motif[,1])
}

## ~ Draw a wrapped multiple sequence alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
wrapAln = function(infile,wrap=80,fex=1,start=1,end=0){
	## Load Data and Setup ##
	aln = read.table(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	seqnum = length(aln[1,])
	seqlen = length(aln[,1])
	if (wrap<1){ wrap = seqlen - start + 1 }
	if (end<start){ end = seqlen }
	blocks = as.integer((end-start)/wrap) + 1
	height = (blocks * (seqnum + 2)) - 1
	## Create empty plot ##
	old_mar = par()$mar
	par(mar=c(0,5,1.1,0))
	#!# png(filename = outfile, width=20*wrap, height=height*30, units = "px")
	par(mar=c(0,5,1.1,0))
	plot(c(0,wrap),c(0,height),type="n",axes=FALSE,ann=FALSE)
	## Plot alignment ##
	for (b in 1:blocks){
		h = height - ((b-1)*(seqnum+2)) - 1
		for (a in 1:seqnum){
			name = colnames(aln)[a]
			text(-0.5,h-a+0.5,colnames(aln)[a],adj=c(1,0.5),xpd=TRUE,font=1,family="sans",cex=fex*min(1,12/seqnum))
			for (i in 1:wrap){
				p = (b-1)*wrap + i + start - 1
				if (p > seqlen){ break }
				if (aln[p,a] == 0 && substr(name,nchar(name)-3,nchar(name)) == ".ppi"){ aa = "-" }
				else { aa = as.character(aln[p,a]) }
				rect(i-1,h-a,i,h-a+1,col=acol[[aa]],border=NA)
				if (aa == "-"){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",col="white",cex=fex)
				}
				if (aa != "-" && nchar(aa) <= 1){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex)
				}
				if (aa != "-" && nchar(aa) > 1){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",srt=90,cex=fex*2/nchar(aa))
				}
				if (i == 1 | i == wrap | as.integer(p/10) == p/10 | p == seqlen){ text(i-0.5,h+0.7,p,adj=c(0.5,1),xpd=TRUE,family="sans",cex=fex*0.8) }
			}
		}
	}
	#!# dev.off()
	par(mar=old_mar)
	height
}
#x#alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=80)

##### ~~~ New wrapped Alignment function ~~~~~~~~~~~~~~~~~~~~~~~~~~ #####
wrapAlnPNG = function(basefile,wrap=80,mode="png"){
	if (mode=="test"){
		wrap=120
		basefile = "SLiMJIM/html/visualisations/PPP1CA.TAI12_HUMAN__Q9H175"
	}
	infile = paste(basefile,".aln.tdt",sep="")
	outfile = paste(basefile,".png",sep="")
	if (mode=="test"){ outfile = "test.png"; mode="png" }
	### ~ [1] ~ Load Data and Setup ~~~~~~~~~~~~~~~~~~~~~~~~~ ###
	aln = read.table(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	## ~ [1a] ~ Alignment stats ~~~~~~~~~~~~~~~~~~ ##
	seqnum = length(aln[1,])				# Number of sequences
	seqlen = length(aln[,1])				# Length of sequences
	namelen = 0							# Length of longest name
	for (a in 1:seqnum){ namelen = max(strwidth(colnames(aln)[a],units="inches"),namelen) }
	if (wrap == 0){ wrap = seqlen }			# No wrap = full length
	blocks = as.integer((seqlen-1)/wrap) + 1		# No of wrapped alignment blocks
	height = (blocks * (seqnum + 2)) - 1		# Height (alignment lines) of wrapped alignment
	## ~ [1b] ~ Calculate plot characteristics ~~~ ##
	fex = 0.8 #0.16 / par()$csi		# Should be 0.8 normally?
	xin = wrap / 8				# Width of plot (inches)
	yin = height / 5				# Height of plot (inches)
	namelen = fex * (namelen + 0.0625)	# Length needed for sequence names (inches)
	if (mode=="size"){ return(c(xin,yin)) }
	## Create empty plot ##
	old_par = par()
	#par(omi=c(0,0,0,0))
	par(mar=c(0.1,5,0.1,0.1))
	if (mode=="png"){ png(filename=outfile, width=xin, height=yin, units="in", res=150) }
	if (mode=="test"){ 
		png(filename=outfile, width=xin, height=yin, units="in", res=150) 
	}
	par(mar=c(0.1,5,0.1,0.1))
	plot(c(0,wrap),c(0,height),type="n",axes=FALSE,ann=FALSE,asp=yin/xin)	#pin=c(xin,yin))
	## Plot alignment ##
	for (b in 1:blocks){
		h = height - ((b-1)*(seqnum+2)) - 1
		for (a in 1:seqnum){
			name = colnames(aln)[a]
			text(-0.5,h-a+0.5,name,adj=c(1,0.5),xpd=TRUE,font=1,family="sans",cex=fex)	#min(1,12/seqnum))
			for (i in 1:wrap){
				p = (b-1)*wrap + i
				if (p > seqlen){ break }
				if (aln[p,a] == 0 && substr(name,nchar(name)-3,nchar(name)) == ".ppi"){ aa = "-" }
				else { aa = as.character(aln[p,a]) }
				rect(i-1,h-a,i,h-a+1,col=acol[[aa]],border=NA)
				if (aa == "-"){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",col="white",cex=fex)
				}
				if (aa != "-" && nchar(aa) <= 1){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex)
				}
				if (aa != "-" && nchar(aa) > 1){
					text(i-0.5,h-a+0.5,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",srt=90,cex=2*fex/nchar(aa))
				}
				if (i == 1 | i == wrap | as.integer(p/10) == p/10 | p == seqlen){ text(i-0.5,h+0.5,p,adj=c(0.5,0.5),xpd=TRUE,family="sans",cex=fex) }
			}
		}
	}
	if (mode=="png"){ 
		dev.off() 
	}
	par(mar=old_par$mar,omi=old_par$omi)
	height
}

## ~ Calculate X and Y coordinates for network nodes, given angle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
getSpokeXY = function(i,inum,offset=0.1){
	### Calculate X and Y from Angle ###
	spangle = 360 / inum
	A1 = spangle * (i-1)
	A = A1
	while (A >= 90){ A = A-90 }
	A = (A / 360) * (2 * pi)	# Convert to rads
	o = sin(A) * (1-offset)
	a = cos(A) * (1-offset)
	if (A1 < 90){ 
		x = o
		y = a
	}
	if (A1 >= 90 & A1 < 180){
		x = a
		y = -o
	}
	if (A1 >= 180 & A1 < 270){
		x = -o
		y = -a
	}
	if (A1 >= 270){
		x = -a
		y = o
	}
	c(x,y)
}

## ~ Plot partial PPI network ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
networkRing = function(infile,h=0.1,w=0.15,fex=1){
	#x#infile = paste(basefile,".ppi.tdt",sep="")
	## Read file and Setup ##
	ppi = read.table(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	spokenum = length(colnames(ppi))
	scale = min(1,10/spokenum)
	h = h * scale
	w = w * scale
	## Empty Plot ##
	old_mar = par()$mar
	par(mar=c(0.2,0.2,0.2,0.2))
	plot(c(-1,1),c(-1,1),type="n",axes=FALSE,ann=FALSE)
	## Draw connections ##
	for (i in 1:spokenum){
		coord = getSpokeXY(i,spokenum,w)
		x = coord[1]
		y = coord[2]
		### Draw connections ###
		for (j in 1:spokenum){
			if (ppi[i,j] == 0){ next }
			coord2 = getSpokeXY(j,spokenum,w)
			lines(c(x,coord2[1]),c(y,coord2[2]),lwd=2)
		}
	}
	### Draw Spokes ###
	for (i in 1:spokenum){
		coord = getSpokeXY(i,spokenum,w)
		x = coord[1]
		y = coord[2]
		if (length(rownames(ppi)) > spokenum){
			rect(x-w,y-h,x+w,y+h,xpd=TRUE,col=ppcol[ppi[spokenum+1,i]+1],lwd=2)
			if (ppi[spokenum+1,i] == 0){
				text(x,y,colnames(ppi)[i],adj=c(0.5,0.5),xpd=TRUE,font=2,cex=fex*1.2*scale,col=soton$col[12])
			}else{
				text(x,y,colnames(ppi)[i],adj=c(0.5,0.5),xpd=TRUE,font=2,cex=fex*1.2*scale)	#,col=soton$col[17])
			}
		}else{
			rect(x-w,y-h,x+w,y+h,xpd=TRUE,col="#979E45",lwd=2)
			text(x,y,colnames(ppi)[i],adj=c(0.5,0.5),xpd=TRUE,font=2,cex=fex*1.2*scale,col=soton$col[1])
		}
	}
	par(mar=old_mar)
}
#x#networkRing(paste(basefile,".ppi.tdt",sep=""))

## ~ Plot partial PPI network ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
nestedNetworkRing = function(infile,h=0.09,w=0.2,fex=1){
	#infile = paste(basefile,".nested.tdt",sep="")
	#h=0.09
	#w=0.2
	## Read file and Setup ##
	ppi = read.table(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	(hubnum = length(colnames(ppi)) - 1)
	(spokenum = length(ppi[,1]))
	iscale = min(1,20/spokenum)		# Size of Outer Ring
	jscale = min(1,10/hubnum)		# Size of Inner Ring
	ih = h * iscale 
	iw = w * iscale
	jh = h * jscale
	jw = w * jscale
	## Empty Plot ##
	old_mar = par()$mar
	par(mar=c(0.2,0.2,0.2,0.2))
	plot(c(-2,2),c(-2,2),type="n",axes=FALSE,ann=FALSE)
	## Draw connections ##
	for (i in 1:spokenum){
		coord = getSpokeXY(i,spokenum,iw)
		x = coord[1] * 2
		y = coord[2] * 2
		### Draw connections ###
		for (j in 1:hubnum){
			if (ppi[i,j+1] == 0){ next }
			coord2 = getSpokeXY(j,hubnum,jw)
			lines(c(x,coord2[1]),c(y,coord2[2]),lwd=2)
		}
	}
	### Draw Hubs ###
	for (j in 1:hubnum){
		coord = getSpokeXY(j,hubnum,jw)
		x = coord[1]
		y = coord[2]
		rect(x-jw,y-jh,x+jw,y+jh,xpd=TRUE,col=soton$col[5],lwd=2)
		text(x,y,colnames(ppi)[j+1],adj=c(0.5,0.5),xpd=TRUE,font=2,cex=fex*1.2*jscale,col=soton$col[1])
	}
	### Draw Spokes ###
	for (i in 1:spokenum){
		coord = getSpokeXY(i,spokenum,iw)
		x = coord[1] * 2
		y = coord[2] * 2
		rect(x-iw,y-ih,x+iw,y+ih,xpd=TRUE,col=soton$col[1],lwd=2)
		text(x,y,ppi[i,1],adj=c(0.5,0.5),xpd=TRUE,font=2,cex=fex*1.2*iscale,col=soton$col[5])
	}
	par(mar=old_mar)
}
#x#networkRing(paste(basefile,".ppi.tdt",sep=""))



## ~ Generation of tree plot from data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
makeTree = function(infile,treetitle="",fex=1){
	### ~ SETUP ~ ###
	tree = read.csv(infile)
	ynum = length(tree$xpos)
	x = 1 / (2*max(tree$xpos))
	y = 1 / (ynum+1)
	tree$fcol = soton$col[1]	#"black"

	### ~ Setup Plot ~ ###
	#!# pngwidth = 1600
	#!# if (length(args) > 3) { pngwidth = as.integer(args[4]); }
	#!# png(filename = outfile, width = pngwidth, height = (ynum+2)*50, units = "px", pointsize=25)
	old_mar = par()$mar
	par(mar=c(0.2,0.2,10,0.2))
	plot(0:1,0:1,type="n",axes=FALSE,ann=FALSE,mar=c(0,1,4,1))
	if (length(treetitle) > 2) { title(main=treetitle); }

	### ~ Draw Tree ~ ###
	for (i in 1:ynum){
		data = tree[i,]
		lines(c(data$xpos*x,data$ancx*x),c(1-data$ypos*y,1-data$ypos*y),col=data$fcol,lwd=2)
		lines(c(data$ancx*x,data$ancx*x),c(1-data$ypos*y,1-data$ancy*y),col=data$fcol,lwd=2)
	}
	### ~ Add Text ~ ###
	for (i in 1:ynum){
		data = tree[i,]
		if (data$nodenum <= ((ynum+1) / 2)){
			text((data$xpos)*x+0.01,1-data$ypos*y,data$name,adj=c(0,0.5),cex=fex*min(1.2,38/ynum),col=data$fcol)
		}else{
			text((data$xpos)*x-0.01,1+(0.45/ynum)-data$ypos*y,data$boot,adj=c(0.8,0),cex=fex*min(1,38/ynum),font=3,col=soton$col[5])
		}
	}
	### ~ Add scale ~ ###
	lines(c(0,0.1*x),c(0,0),col=soton$col[7],lwd=2)
	lines(c(0,0),c(0,-0.005),col=soton$col[7],lwd=2)
	lines(c(0.1*x,0.1*x),c(0,-0.005),col=soton$col[7],lwd=2)
	text(0,-0.01,"0",adj=c(0.5,1),cex=fex*min(0.8,38/ynum),col=soton$col[7],xpd=NA)
	text(0.1*x,-0.01,"0.1",adj=c(0.5,1),cex=fex*min(0.8,38/ynum),col=soton$col[7],xpd=NA)
	par(mar=old_mar)
	#!# dev.off()
}

## ~ Generation of heatmap from data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
heatMap = function(infile,fex=1){
	## Read file and Setup ##
	gb = 255 - 144
	dist = read.table(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	names = dist[,1]
	n = 1 / length(names)
	## Empty Plot ##
	old_mar = par()$mar
	par(mar=c(2,0.2,8,0.2))
	plot(0:1,0:1,type="n",axes=FALSE,ann=FALSE)
	## Squares ##
	for (i in 1:length(names)){
		#text(-0.01,1-(i-0.5)*n,dist[i,1],adj=c(1,0.5),xpd=TRUE,font=1,family="sans",cex=0.7)
		text((i-0.5)*n,1.01,dist[i,1],adj=c(0,0.5),xpd=TRUE,font=1,family="sans",cex=fex,srt=90)
		for (j in 1:length(names)){
			D = 1 - dist[i,j+1]
			C = (255 - (D * gb)) / 255
			rect((i-1)*n,1-(j-1)*n,i*n,1-j*n,col=rgb(1-D,C,C),xpd=TRUE,lwd=2)
			text((i-0.5)*n,1-(j-0.5)*n,as.integer(100*(D+0.005)),adj=c(0.5,0.5),xpd=TRUE,cex=fex*min(1,12/length(names)))
		}
	}
	par(mar=old_mar)
}

############# ::: PARAMETER SETUP ::: ###################################################################################

## ~ Read commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
slimjim = "arse"
args <- commandArgs(TRUE)
if (length(args) > 1){
	basefile = args[1]
	slimjim = args[2]
}

############# ::: MAIN RUN - GRAPHICS PLOT ::: ##########################################################################

########## ~~~ MOTIF MULTIPANEL PLOT ~~ ##############
multiPanelMotif = function(basefile){
	## ~ Setup multipanel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	#alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=80)
	png(filename=paste(basefile,".png",sep=""), width=2400, height=1600, units = "px", pointsize=12)
	#png(filename=paste(basefile,".png",sep=""), width=2400, height=1600, units = "px", pointsize=32)
	#panels = c(1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,3,3,3,4,4,4,5,5,5,3,3,3,4,4,4,5,5,5)
	panels = c(1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,2,2,2,3,4,4,4,5,5,5,5,5,3,4,4,4,5,5,5,5,5,3,4,4,4,5,5,5,5,5)
	layout(matrix(panels,byrow=TRUE,nrow=5))
	#layout.show(11)
	fex = 2.67

	## ~ Plot alignment surrounding motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=0,fex=fex)
	
	## ~ Plot sequence logo of motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	mlen = seqLogo(paste(basefile,".profile.tdt",sep=""),fex=fex)
	
	## ~ Plot tree for motif-containing spokes (based on GABLAM) ~~~~~~~~~ ##
	makeTree(paste(basefile,".tree.csv",sep=""),fex=fex)
	
	## ~ Make Heatmap (based on GABLAM) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	heatMap(paste(basefile,".heatmap.tdt",sep=""),fex=fex)
	
	## ~ Plot PPI network for motif-containing spokes ~~~~~~~~~~~~~~~~~~~~ ##
	networkRing(paste(basefile,".ppi.tdt",sep=""),fex=fex)
	dev.off()
}
if (slimjim == "motif"){ multiPanelMotif(basefile) }

########## ~~~ PPI MULTIPANEL PLOT ~~ ##############
multiPanelPPI = function(basefile){
	## ~ Setup multipanel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	png(filename=paste(basefile,".png",sep=""), width=1500, height=2400, units = "px", pointsize=12)
	panels = c(1,1,1,1,1,2,3,3,3,3)
	layout(matrix(panels,byrow=TRUE,nrow=2))
	#layout.show(11)
	## ~ Plot PPI network for motif-containing spokes ~~~~~~~~~~~~~~~~~~~~ ##
	networkRing(paste(basefile,".ppi.tdt",sep=""),fex=3)
	## ~ Plot tree for motif-containing spokes (based on GABLAM) ~~~~~~~~~ ##
	makeTree(paste(basefile,".tree.csv",sep=""),fex=3)
	## ~ Make Heatmap (based on GABLAM) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	heatMap(paste(basefile,".heatmap.tdt",sep=""),fex=2)
	dev.off()
}
if (slimjim == "interactome"){ multiPanelPPI(basefile) }


######### MISCELLANEOUS EXTRA CODE #####################
#plot(c(0,20),c(0,1))
#for (i in 1:21){ rect(i-1,0,i,1,col=soton$col[i]) }
#for (i in 1:21){ text(i-0.5,0.5,i) }

if (slimjim == "test"){
	basefile = "SLiMJIM\\html\\visualisations\\PPP1CA.1.KR-1-IV-1-F-0-AS"
	multiPanelMotif(basefile)
	basefile = "SLiMJIM\\html\\visualisations\\PPP1CA.2.S-0-D-1-DE-1-DE"
	multiPanelMotif(basefile)
	basefile = "SLiMJIM\\html\\visualisations\\PPP1CA.3.KR-0-R-0-V-1-F-0-A"
	multiPanelMotif(basefile)
	basefile = "SLiMJIM\\html\\visualisations\\PPP1CA.4.KR-0-V-1-F-0-AS"
	multiPanelMotif(basefile)
	basefile = "SLiMJIM\\html\\visualisations\\PPP1CA.5.KR-0-KR-0-V-1-F"
	multiPanelMotif(basefile)
	basefile = "SLiMJIM\\html\\visualisations\\PPP1CA.6.KR-1-FV-1-F-1-DE"
	multiPanelMotif(basefile)
	basefile = "SLiMJIM\\html\\visualisations\\PPP1CA.7.V-1-F-1-D-0-E"
	multiPanelMotif(basefile)
	basefile = "SLiMJIM\\html\\visualisations\\PPP1CA.interactome"
	multiPanelPPI(basefile)

	basefile = "SLiMJIM\\html\\visualisations\\PPP1CA.YLPM1_HUMAN__P49750"
	alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=80)
	png(filename=paste(basefile,".png",sep=""), width=1800, height=20*alnh, units = "px", pointsize=12)
	alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=80,fex=2.67)
	dev.off()
}


#wrapAlnPNG("SLiMJIM/html/visualisations/PPP1CA.TAI12_HUMAN__Q9H175",0,"size")
#wrapAlnPNG(mode="test")

######## ~~~ RELATIVE CONSERVATION AND DISORDER PLOT ~~~~ #########################
relConsPlot = function(infile,start=1,end=0,fex=1){
	## Load Data and Setup ##
	rel = read.table(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE,as.is=TRUE)
	seqlen = length(rel[,1])
	if (end<1){ end = seqlen }
	wrap = end - start
	## Create empty plot ##
	old_mar = par()$mar
	par(mar=c(0,5,1.1,0))
	#!# png(filename = outfile, width=20*wrap, height=height*30, units = "px")
	par(mar=c(0,5,1.1,0))
	plot(c(0,wrap),c(-1.5,1.5),type="n",axes=FALSE,ann=FALSE)
	axis(side=2,lwd=2,cex=fex)
	## Plot alignment ##
	for (i in 0:wrap){
		r = start + i
		g = 0.1	# Size of gap each side of bar
		dis = rel[r,"IUPred"]
		ord = rel[r,"Disorder"]
		pos = rel[r,"Pos"]
		aa = rel[r,2]
		cons = rel[r,"RelVNE"]
		# Plot relcons, then disorder, then sequence on top. #
		if (cons<0){
			rect(i+g,0,i+1-g,cons,lwd=2,border=soton$col[5],col=soton$col[9])
		}else{
			rect(i+g,0,i+1-g,cons,lwd=2,border=soton$col[5],col=soton$col[5])
		}
		if (ord == "Order"){
			rect(i+g,0,i+1-g,dis,lwd=2,border=soton$col[1],col=soton$col[2])
		}else{
			rect(i+g,0,i+1-g,dis,lwd=2,border=soton$col[1])
		}
		if (cons>0 & cons<dis & ord=="Order"){ rect(i+g,0,i+1-g,cons,lwd=2,border=soton$col[5]) }
		# Plot sequence and position #
		rect(i,-0.01,i+1,-0.11,col="white",border="white")
		rect(i,-0.015,i+1,-0.105,col=acol[[aa]],border=NA)
		text(i+0.5,-0.06,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex)
		if (i == 0 | i == wrap | as.integer(r/10) == r/10 | r == seqlen){ text(i+0.5,-0.16,r,adj=c(0.5,0.5),xpd=TRUE,family="sans",cex=fex) }
	}
	par(mar=old_mar)
}

######## ~~~ RELATIVE CONSERVATION, DISORDER AND ALIGNMENT PLOT ~~~~ #########################
relConsAln = function(basefile,start=1,end=0,fex=1,lex=2){
	## Load Data and Setup ##
	rel = read.table(paste(basefile,".rel.tdt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE,as.is=TRUE)
	aln = read.table(paste(basefile,".aln.tdt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE,as.is=TRUE)
	seqlen = length(rel[,1])
	if (end<1){ end = seqlen }
	wrap = end - start
	## Create empty plot ##
	old_mar = par()$mar
	par(mar=c(0,5,1.1,0),ljoin="mitre")
	#!# png(filename = outfile, width=20*wrap, height=height*30, units = "px")
	par(mar=c(0,5,1.1,0))
	plot(c(0,wrap),c(-1.5,1.5),type="n",axes=FALSE,ann=FALSE)
	g = 0.1	# Size of gap each side of bar
	axis(side=2,lwd=2,cex=fex,pos=-(2*g))
        rect(-(4*g)-1,-1.6,-(4*g),1.6,col="white",xpd=TRUE,border=NA)
 	for (y in c(1.5,1.0,0.5,0.0,-0.5,-1.0,-1.5)){
 	        text(-(5*g),y,y,adj=c(1,0.5),xpd=TRUE,family="sans",cex=0.8*fex)
        }
	## Plot alignment ##
	for (i in 0:wrap){
		r = start + i
		if(r>seqlen){ break }
		dis = rel[r,"IUPred"]
		ord = rel[r,"Disorder"]
		pos = rel[r,"Pos"]
		aa = rel[r,2]
		cons = rel[r,"RelVNE"]
		# Plot relcons, then disorder, then sequence on top. #
		if (cons<0){
			rect(i+g,0,i+1-g,cons,lwd=lex,border=soton$col[5],col=soton$col[9])
		}else{
			rect(i+g,0,i+1-g,cons,lwd=lex,border=soton$col[5],col=soton$col[5])
		}
		if (ord == "Order"){
			rect(i+g,0,i+1-g,dis,lwd=lex,border=soton$col[1],col=soton$col[2])
		}else{
			rect(i+g,0,i+1-g,dis,lwd=lex,border=soton$col[1])
		}
		if (cons>0 & cons<dis & ord=="Order"){ rect(i+g,0,i+1-g,cons,lwd=lex,border=soton$col[5]) }
		# Plot sequence and position #
		rect(i,-0.01,i+1,-0.11,col="white",border="white")
		rect(i,-0.015,i+1,-0.105,col=acol[[aa]],border=NA)
		text(i+0.5,-0.06,aa,adj=c(0.5,0.5),xpd=FALSE,font=2,family="mono",cex=fex)
		if (i == 0 | i == wrap | as.integer(r/10) == r/10 | r == seqlen){ text(i+0.5,-0.16,r,adj=c(0.5,0.5),xpd=TRUE,family="sans",cex=fex) }
		# Add Motif #
		if (aln[r,1] != "-" & as.character(aln[r,1]) != "0"){
            	h = min(1.5,max(cons,dis) + 0.05)
			rect(i,h,i+1,h+0.09,col=acol[[aln[r,1]]],border=NA,xpd=TRUE)
		   	text(i+0.5,h+0.045,aln[r,1],adj=c(0.5,0.5),xpd=TRUE,font=2,family="mono",cex=fex)
            }

	}
	## Add key ##
	k = 0.02
	rect(k,-1.0-k,11+k,-1.125+k,lwd=2,border=NA,col="white")
	rect(k,-1.0-k,1+k,-1.125+k,lwd=2,border=soton$col[5],col=soton$col[5])
	text(1.5+k,-1.0625,"Conserved",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	rect(k,-1.125-k,11+k,-1.25+k,lwd=2,border=NA,col="white")
	rect(k,-1.125-k,1+k,-1.25+k,lwd=2,border=soton$col[5],col=soton$col[9])
	text(1.5+k,-1.1875,"Unconserved",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	rect(k,-1.25-k,11+k,-1.375+k,lwd=2,border=NA,col="white")
	rect(k,-1.25-k,1+k,-1.375+k,lwd=2,border=soton$col[1],col=soton$col[2])
	text(1.5+k,-1.3125,"Ordered",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	rect(k,-1.375-k,11+k,-1.5+k,lwd=2,border=NA,col="white")
	rect(k,-1.375-k,1+k,-1.5+k,lwd=2,border=soton$col[1])
	text(1.5+k,-1.4375,"Disordered",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	par(mar=old_mar)
}

########## ~~~ SPOKE ALIGN PLOT ~~ ##############

testSpokeAln = function(basefile){
	## ~ Single panel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	aln = read.table(paste(basefile,".aln.tdt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE)
	seqnum = length(aln[1,])
	seqlen = length(aln[,1])
	wrap = 80
	if (wrap<1){ wrap = seqlen }
	blocks = as.integer((seqlen-1)/wrap) + 1
	alnh = (blocks * (seqnum + 2)) - 1
	png(filename=paste(basefile,".png",sep=""), width=20+20*wrap, height=20*alnh, units = "px", pointsize=12)
	alnh = wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=wrap,fex=1.5)
	dev.off()
	png(filename=paste(basefile,".rel.png",sep=""), width=20+20*seqlen, height=600, units = "px", pointsize=12)
	relConsAln(basefile,fex=0.8,lex=2)
	dev.off()


	## ~ Setup multipanel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	wrap = 80
	if (wrap<1){ wrap = seqlen }
	blocks = as.integer((seqlen-1)/wrap) + 1
	alnh = (blocks * (seqnum + 2)) - 1
	png(filename=paste(basefile,".2.png",sep=""), width=20+20*wrap, height=20*alnh+400*blocks, units = "px", pointsize=12)
	panels = c(rep(1,blocks),2:(blocks+1))
	layout(matrix(panels,byrow=TRUE,ncol=1))
	#par(mfcol=c(blocks*2,1))
	wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=wrap,fex=1.5)
	for (b in 1:blocks){
		relConsAln(basefile,start=1+(b-1)*wrap,end=b*wrap)	#!# Add splitting and wrapping #!#
	}
	dev.off()

	## ~ Setup multipanel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	wrap = 0
	if (wrap<1){ wrap = seqlen }
	blocks = as.integer((seqlen-1)/wrap) + 1
	alnh = (blocks * (seqnum + 2)) - 1
	png(filename=paste(basefile,".3.png",sep=""), width=20+20*wrap, height=20*alnh+400*blocks, units = "px", pointsize=12)
	panels = c(1:2)
	layout(matrix(panels,byrow=TRUE,ncol=1))
	#par(mfcol=c(blocks*2,1))
	wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=wrap,fex=1.5)
	for (b in 1:blocks){
		relConsAln(basefile,start=1+(b-1)*wrap,end=b*wrap)	#!# Add splitting and wrapping #!#
	}
	dev.off()
}
spokeAln = function(basefile){
	## ~ Single panel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	aln = read.table(paste(basefile,".aln.tdt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE)
	seqnum = length(aln[1,])
	seqlen = length(aln[,1])
	start = 1
	end = 110
	wrap = 110
	while (start < seqlen){
		pngfile = paste(basefile,".",preZero(start),"-",preZero(end),".png",sep="")
		alnh = seqnum + 2
		png(filename=pngfile, width=2400, height=min(1200,80*alnh), units = "px", pointsize=12)
		panels = c(rep(1,as.integer((seqnum-1)/25)+1),rep(2,2))
		layout(matrix(panels,byrow=TRUE,ncol=1))
		wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=wrap,fex=3,start=start,end=end)
		relConsAln(basefile,fex=3,lex=2,start=start,end=end)
		dev.off()
		if (end > seqlen){ break }
		start = start + 100
		end = end + 100
	}
}
if (slimjim == "spokealn"){ spokeAln(basefile) }

nestedPNG = function(basefile){
	png(filename=paste(basefile,".png",sep=""), width=2400, height=2400, units = "px", pointsize=12)
	nestedNetworkRing(paste(basefile,".nested.tdt",sep=""),fex=3)
	dev.off()
}
if (slimjim == "nested"){ nestedPNG(basefile) }


