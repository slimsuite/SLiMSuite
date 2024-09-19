########################################################
### RJE_PLOTS: Reusable R plotting functions   ~~~~~ ###
### VERSION: 0.1.0                             ~~~~~ ###
### LAST EDIT: 27/09/23                        ~~~~~ ###
### AUTHORS: Richard Edwards 2023              ~~~~~ ###
### CONTACT: rich.edwards@uwa.edu.au           ~~~~~ ###
########################################################

# This script is for mapping/updating sequence identifiers

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial version for assembly QV analysis
version = "v0.1.0"

####################################### ::: SUMMARY ::: ############################################
#i# vioPlot(plotdb,plotfield) = Violin plot function from DepthKopy, dividing by dataset.

####################################### ::: TODO ::: ############################################
# [ ] : Generalise the plotting function

### ~~~~ Violin plot function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
vioPlot <- function(plotdb,plotfield){
  #i# Violin plots by Method
  #!# Add labels of n=X for each Dataset
  labels = levels(factor(plotdb$Dataset))
  labels = paste("n =",table(factor(plotdb$Dataset)))
  nx = 1:length(labels)
  ny = rep(0,length(labels))
  p <- ggplot(plotdb, aes(factor(Dataset), .data[[plotfield]]), palette = palette) + 
    aes(fill = factor(Dataset)) +
    geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) + 
    xlab("Dataset") + ylab(plotfield) +
    labs(fill="Dataset") +
    #geom_text(x=nx, y=ny, label=labels, vjust=1) +
    theme_bw(base_size = settings$pointsize)
  
  if(settings$ggstatsplot){
    splotdb = plotdb
    splotdb$plotfield = splotdb[[plotfield]]
    p <- ggbetweenstats(
      data = splotdb,
      results.subtitle = length(labels) <= 4 & settings$sigdif,
      centrality.label.args = list(size  = settings$pointsize/4, nudge_y = 0.5, alpha = 0.75),
      pairwise.comparisons = settings$sigdif,
      x = Dataset,
      y = plotfield
    ) + ylim(0,median(splotdb$plotfield)*settings$cnmax) +
      #geom_hline(yintercept=1, color="black") +
      labs(
        #title = "DepthSizer accuracy versus reference assembly size",
        x = "Dataset",
        y = plotfield,
        caption = NULL
      ) +
      theme(text = element_text(size=settings$pointsize))
  }
  
  nlevel = length(levels(factor(plotdb$Dataset)))
  if(nlevel > 8){
    p <- p + scale_colour_manual(values=rep(brewer.pal(12,"Paired"),times=as.integer(nlevel/12)+1)[1:nlevel])
  }
  
  if("Other" %in% levels(factor(plotdb$Dataset))){
    p <- p + theme(axis.text.x = element_text(angle = 90))
  }
  
  if(plotfield %in% c("CN","SelfK")){
    p <- p + geom_hline(yintercept=1, color="steelblue", linetype="dashed")
  }
  if(plotfield %in% c("MeanX","MedX","ModeX","DensX")){
    p <- p + geom_hline(yintercept=scdepth, color="steelblue", linetype="dashed")
  }
  
  return(p)
}
