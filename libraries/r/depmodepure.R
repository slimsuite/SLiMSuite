#i# v0.1 - Modified to use adjust=12 for smoother curve.
# Example use:
# Rscript depmodepure.R tmpfile

#i# Load data from commandline argument
infile = commandArgs(TRUE)[1]
dephist = read.delim(infile, sep="", header=FALSE)
#i# Extract depmethod
colnames(dephist) = c("X")
#i# Calculate "mode" of density distribution
depvec = dephist$X
depdens = density(depvec,n=2048,adjust=12)
depmode = depdens$x[depdens$y == max(depdens$y)]
#i# Write to STDOUT
writeLines(as.character(depmode))

#i# Straight mode of depth vector
getmode <- function(v) {
  uniqv <- unique(v)
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}
#i# Density plotting function
densModePlot <- function(depvec,bdens,plotmain){
  centre = getmode(depvec)
  meanx = mean(depvec)
  dmode=bdens$x[bdens$y == max(bdens$y)]

  plot(bdens,xlim=c(0,dmode*2), lwd = 2,col="blue",
       main=plotmain,ylab="Density",xlab="Score",
       ylim=c(0,max(bdens$y)))

  abline(v=centre,col="grey")
  text(centre-2,0,round(centre,2),adj=c(1,0),col="grey")

  abline(v=dmode,col="steelblue")
  text(dmode+2,0,round(dmode,2),adj=0,col="steelblue")
  abline(h=max(bdens$y),col="steelblue")

  abline(v=meanx,col="red")
  text(meanx+2,max(bdens$y)/2,round(meanx,2),adj=0,col="red")

}

pngfile = paste0(infile,".density.png")
png(pngfile,width=1200,height=900,pointsize=16)
densModePlot(depvec,depdens,paste(infile,"raw density profile"))
dev.off()

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# depmode = float(rje.chomp(os.open('Rscript {0}depmodepure.R {1} {2}'.format(rdir, tmpfile)).readlines()[0]))
