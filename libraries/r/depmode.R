# Example use:
# Rscript depmode.R ../data/2020-09-28.Diploidocus/e_coli.dephist.tdt depmethod

#i# Load data from commandline argument
infile = commandArgs(TRUE)[1]
depmethod = commandArgs(TRUE)[2]
dephist = read.delim(infile, sep="\t")
#i# Extract depmethod
dephist = dephist[dephist$Method == depmethod,]
#i# Calculate "mode" of density distribution
depvec = rep(dephist$X,dephist$n)
depdens = density(depvec)
depmode = depdens$x[depdens$y == max(depdens$y)]
#i# Write to STDOUT
writeLines(as.character(depmode))

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# depmode = float(rje.chomp(os.open('Rscript {0}depmode.R {1} {2}'.format(rdir, depfile, depmethod)).readlines()[0]))
