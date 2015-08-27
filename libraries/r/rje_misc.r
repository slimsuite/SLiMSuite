######################################################
### GENERAL USAGE FUNCTIONS FOR RJE R SCRIPTS ~~~~ ###
### AUTHOR: Dr Richard Edwards 2008 ~~~~~~~~~~~~~~ ###
### CONTACT: r.edwards@soton.ac.uk ~~~~~~~~~~~~~~~ ###
######################################################

############# ::: GENERAL FUNCTIONS ::: ################

### ~ Adding leading zeros to numbers and returning as string ~ ###
preZero = function(n,digits=5){
	nstr = as.character(n)
	while (nchar(nstr)<digits){ nstr = paste("0",nstr,sep="") }
	return(nstr)
}

