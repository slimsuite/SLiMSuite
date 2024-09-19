######################################################
### GO ENRICHMENT PNG SCRIPT TO ACCOMPANY APHID.PY ###
### AUTHOR: Dr Richard Edwards 2008 ~~~~~~~~~~~~~~ ###
### CONTACT: r.edwards@soton.ac.uk ~~~~~~~~~~~~~~~ ###
######################################################

############# ::: GENERAL SETUP ::: ##################
args <- commandArgs(TRUE)
goset = args[1]
basefile = args[2]
pcut = 0.05
if (length(args) > 2){ pcut = as.numeric(args[3]) }
#goset = 'goslim_power'
#goset = 'goslim_generic'


############# ::: MAKE GRAPHIC FUNCTION ::: ##################
go_graph = function(sample,pcut=0.05,fex=1,usebonf=FALSE){
     # Sample = Sample for analysis
     # pcut = Probability cut-off
     # usebonf = Perform Bonferroni correction
     ## ~ [1] ~ Setup ~~~~~ ##
     panels = c(1,4,4,2,5,5,3,6,6)
	 layout(matrix(panels,byrow=TRUE,nrow=3))
     sc = 2	# Colour index
     ## ~ [2] ~ Vs Combined Data ~~~~~~~~ ##
     bonf=1
     comp = "vs Combined Data"
     if (usebonf){ (bonf = length(which(aphid["Pingu"] > 0 & aphid[goset] > 0))) }
     for (gotype in c("bp","mf","cc")) {

		logo = paste("logo_",sample,sep="")
		aphid[logo] = log(aphid[sample]/aphid[paste("e_",sample,sep="")],base=2)

		(prows = which(aphid$Type==gotype & (aphid[paste("p_",sample,sep="")] * bonf) <= pcut))
        if (length(prows)<1){
            plot(c(0,1),c(-2,2),type="n",axes=FALSE,ann=FALSE,main=paste(sample,gotype,comp),cex.axis=fex,cex.names=fex,cex.lab=fex,cex.main=fex,cex.sub=fex)
            next
        }
		#if (length(prows) > 20){
		#	(prows = which(order(aphid[,paste("p_",sample,sep="")])<=20))
		#}
		barplot(aphid[logo][prows,],ylim=c(-2,2),
			col=soton$col[sc],
			border=soton$col[1],
			main=paste(sample,gotype,comp),	#"p<=",pcut),
			ylab="Log2 GO Term enrichment",cex.axis=fex,cex.names=fex,cex.lab=fex,cex.main=fex,cex.sub=fex)
		#text(c(1:length(aphid$Desc[prows])),aphid[logo][prows,],aphid$Desc[prows],srt=90)
		(blabels = as.character(aphid$Desc[prows]))
		(bsig = aphid[paste("p_",sample,sep="")][prows,] * bonf)

		alabels = blabels			# Add asterisks to these labels
		for (i in 1:length(blabels)) {
			alabels[i] = ifelse(nchar(blabels[i])>40,paste(substr(blabels[i],1,40),"...",sep=""),blabels[i])
		}
		(plabels = alabels)
		stars = c('*','**','***')
		scuts = c(0.5,0.01,0.001)
		for (i in 1:3) {
		   (asig = which(bsig <= scuts[i]))
		   plabels[asig] = paste(alabels[asig],stars[i])
		}
		text(c(1:length(aphid$Desc[prows]))*1.2-0.6,c(0),plabels,srt=90,xpd=TRUE)

		sc = sc + 1
	}
	sc = sc + 3
    ## ~ [3] ~ Vs EnsEMBL ~~~~~~~~~~~ ##
    comp = "vs EnsEMBL"
    if (usebonf){ (bonf = length(which(aphid["ense_Activated"] > 0))) }
	for (gotype in c("bp","mf","cc")) {

		logo = paste("logo_",sample,sep="")
		aphid[logo] = log(aphid[sample]/aphid[paste("ense_",sample,sep="")],base=2)

		(prows = which(aphid$Type==gotype & (aphid[paste("ensp_",sample,sep="")] * bonf) <= pcut))
        if (length(prows)<1){
            plot(c(0,1),c(0,6),type="n",axes=FALSE,ann=FALSE,main=paste(sample,gotype,comp),cex.axis=fex,cex.names=fex,cex.lab=fex,cex.main=fex,cex.sub=fex)
            next
        }

		#if (length(prows) > 20){
		#	(prows = which(order(aphid[,paste("ensp_",sample,sep="")])<=20))
		#}
		barplot(aphid[logo][prows,],ylim=c(0,6),
			col=soton$col[sc],
			border=soton$col[1],
			main=paste(sample,gotype,comp),	#"p<=",pcut),
			ylab="Log2 GO Term enrichment",cex.axis=fex,cex.names=fex,cex.lab=fex,cex.main=fex,cex.sub=fex)
		#text(c(1:length(aphid$Desc[prows])),aphid[logo][prows,],aphid$Desc[prows],srt=90)
		(blabels = as.character(aphid$Desc[prows]))
		(bsig = aphid[paste("ensp_",sample,sep="")][prows,] * bonf)

		alabels = blabels			# Add asterisks to these labels
		for (i in 1:length(blabels)) {
			alabels[i] = ifelse(nchar(blabels[i])>40,paste(substr(blabels[i],1,40),"...",sep=""),blabels[i])
		}
		(plabels = alabels)
		stars = c('*','**','***')
		scuts = c(0.5,0.01,0.001)
		for (i in 1:3) {
		   (asig = which(bsig <= scuts[i]))
		   plabels[asig] = paste(alabels[asig],stars[i])
		}
		text(c(1:length(aphid$Desc[prows]))*1.2-0.6,c(3),plabels,srt=90,xpd=TRUE,cex=fex)

		sc = sc + 1
	}
	sc = sc + 3
}###END###



############### ::: READ IN DATA ::: ################
aphid = read.table(paste(basefile,goset,"tdt",sep="."),sep="\t",header=TRUE)
dim(aphid)
rownames(aphid) = as.character(aphid$GO)

colnum = length(colnames(aphid))
samplex = as.integer((colnum - 7)/5)
samples = colnames(aphid)[8:(7+samplex)]



############## ::: MAKE FIGURES ::: ################
for (sample in samples){
	png(filename=paste(basefile,sample,goset,"png",sep="."), width=2400, height=1500, units = "px", pointsize=12)
    go_graph(sample,pcut,fex=2.5)
	dev.off()
}
