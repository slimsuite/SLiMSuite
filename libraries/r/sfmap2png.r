########################################################
### COMPOSITE IMAGE SCRIPT TO ACCOMPANY SFMAP2PNG.PY ###
### AUTHOR: Dr Richard Edwards 2008	~~~~~~~~~~~~~~~~ ###
### CONTACT: r.edwards@soton.ac.uk ~~~~~~~~~~~~~~~~~ ###
########################################################

############# ::: GENERAL SETUP ::: #####################################################################################

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
	par(mar=c(0,10,1.1,0))
	#!# png(filename = outfile, width=20*wrap, height=height*30, units = "px")
	par(mar=c(0,10,1.1,0))
	plot(c(0,wrap-1),c(0,height),type="n",axes=FALSE,ann=FALSE)
	## Plot alignment ##
	for (b in 1:blocks){
		h = height - ((b-1)*(seqnum+2)) - 1
		for (a in 1:seqnum){
			name = colnames(aln)[a]
			text(-0.5,h-a+0.5,colnames(aln)[a],adj=c(1,0.5),xpd=TRUE,font=1,family="sans",cex=fex*min(0.9,12/seqnum))
			for (i in 1:wrap){
				p = (b-1)*wrap + i + start - 1
				if (p > seqlen){ break }
				if (aln[p,a] == 0 && a == 1){ aa = "-" }
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

############# ::: PARAMETER SETUP ::: ###################################################################################

## ~ Read commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
args <- commandArgs(TRUE)
basefile = args[2]

############# ::: MAIN RUN - GRAPHICS PLOT ::: ##########################################################################
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
	par(mar=c(0,10,1.1,0),ljoin="mitre")
	#!# png(filename = outfile, width=20*wrap, height=height*30, units = "px")
	par(mar=c(0,10,1.1,0))
	plot(c(0,wrap),c(-1.5,1.5),type="n",axes=FALSE,ann=FALSE)
	g = 0.1	# Size of gap each side of bar
	axis(side=2,lwd=2,cex=fex,pos=-(2*g))
        rect(-(4*g)-1,-1.6,-(4*g),1.6,col="white",xpd=TRUE,border=NA)
 	for (y in c(1.5,1.0,0.5,0.0,-0.5,-1.0,-1.5)){
 	        text(-(5*g),y,y,adj=c(1,0.5),xpd=TRUE,family="sans",cex=fex)
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
		# Add Motif #
                if (aln[r,1] != "-" & as.character(aln[r,1]) != "0"){
  	                h = min(1.5,max(cons,dis) + 0.05)
			rect(i,h,i+1,h+0.09,col=acol[[aln[r,1]]],border=NA,xpd=TRUE)
		   	text(i+0.5,h+0.045,aln[r,1],adj=c(0.5,0.5),xpd=TRUE,font=2,family="mono",cex=fex)
                }

	}
	for (i in 0:wrap){
		r = start + i
		if(r>seqlen){ break }
		if (i == 0 | i == wrap | as.integer(r/10) == r/10){ text(i+0.5,-0.16,r,adj=c(0.5,0.5),xpd=TRUE,family="sans",cex=fex) }
	}
	## Add key ##
	k = 0.02
	#rect(k,-1.0+k,11+k,-1.125-k,lwd=2,border=NA,col="white")
	#rect(k,-1.125+k,11+k,-1.25-k,lwd=2,border=NA,col="white")
	#rect(k,-1.25+k,11+k,-1.375-k,lwd=2,border=NA,col="white")
	#rect(k,-1.375+k,11+k,-1.5-k,lwd=2,border=NA,col="white")
	rect(0.5,-1.0+k,11.5,-1.5-k,lwd=2,border=soton$col[6],col="white")
	rect(1,-1.0-k,2,-1.125+k,lwd=2,border=soton$col[5],col=soton$col[5])
	text(2.5,-1.0625,"Conserved",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	rect(1,-1.125-k,2,-1.25+k,lwd=2,border=soton$col[5],col=soton$col[9])
	text(2.5,-1.1875,"Unconserved",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	rect(1,-1.25-k,2,-1.375+k,lwd=2,border=soton$col[1],col=soton$col[2])
	text(2.5,-1.3125,"Ordered",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	rect(1,-1.375-k,2,-1.5+k,lwd=2,border=soton$col[1])
	text(2.5,-1.4375,"Disordered",adj=c(0,0.5),xpd=TRUE,family="sans",font=2,cex=fex)
	par(mar=old_mar)
}

########## ~~~ SPOKE ALIGN PLOT ~~ ##############
spokeAln = function(basefile,method="png"){
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
		rel = 7
		panels = c(rep(1,as.integer((alnh-1)/4)+1),rep(2,rel))
		ph = as.integer(40.0 * alnh * length(panels) / (length(panels)-rel))
		#png(filename=pngfile, width=2400, height=min(1200,80*alnh), units = "px", pointsize=12)
		if (method == "png"){
		   png(filename=pngfile, width=2400, height=ph, units = "px", pointsize=12)
                }else{
		      CairoPNG(filename=pngfile, width=2400, height=ph, pointsize=12)
                }
		#panels = c(rep(1,as.integer((alnh-1)/8)+1),rep(2,4))
		layout(matrix(panels,byrow=TRUE,ncol=1))
		wrapAln(paste(basefile,".aln.tdt",sep=""),wrap=wrap,fex=2.5,start=start,end=end)
		relConsAln(basefile,fex=2.5,lex=2,start=start,end=end)
		dev.off()
		if (end > seqlen){ break }
		start = start + 100
		end = end + 100
	}
}
spokeAln(basefile)

