########################################################
### SynBad Plot Functions                      ~~~~~ ###
### VERSION: 0.0.0                             ~~~~~ ###
### LAST EDIT: 03/01/21                        ~~~~~ ###
### AUTHORS: Richard Edwards 2016              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for genomic plotting functions to be used by SynBad analysis Rmd.
# In future, a synbad.R script might generate some graphical outputs automatically.

################# ::: HISTORY ::: ######################
# v0.0.0 : Initial version based on SynBad 2021-11-24 analysis.

################# ::: TO DO ::: ######################
# [ ] : Transfer and generalise functions.

################# ::: SETUP ::: ######################

if(! "settings" %in% ls()){
  settings <- list(outlog=stdout())
}

if(! "logWrite" %in% ls()){
    if(! outlog %in% names(settings)){ settings$outlog <- stdout() }
    logWrite <- function(logstr){
      writeLines(paste0("[",date(),"] ",logstr),con=settings$outlog)
    }
}

################# ::: SYNBAD FUNCTIONS ::: ######################

### ~ Function template ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# > D : data frame with Prefix, SeqName, Start, End, and Strand
# > logwrite=TRUE : Whether to summarise to logWrite()
newfunc <- function(D,logwrite=TRUE){
    if(nrow(D) < 1)){
       abort("No data!")
    }
    # Summarise
    if(logwrite & "logWrite" %in% ls()){
      logWrite(paste(nrow(D),"rows"))
    }
    return(D)
}


