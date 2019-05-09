#i# See the `main.R` file for version, contact and license information.
#i# The `main.R` script loads libraries and contains the initial parameter settings and functions.
#i# This script is loaded by `main.R` and functions as a convenient place for the primary SNP filtering function.
#i# This in primarily to make it easier to find and edit the plotting code as new functions are added.

## ~ FreqTable SNP filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
filterSNPs = function(appdata,freqdb){
  snpdata = freqdb
  ### Subsetting data
  if(appdata$freqfilter == "Fixed"){
    snpdata = snpdata[snpdata$MajFreq >= (1.0-appdata$fixfreq),]  
  }
  if(appdata$freqfilter == "Lost"){
    snpdata = snpdata[snpdata$MajFreq <= appdata$fixfreq,]  
  }
  if(appdata$freqfilter == "Polymorphic"){
    snpdata = snpdata[snpdata$MajFreq > appdata$fixfreq & snpdata$MajFreq < (1.0-appdata$fixfreq),]  
  }
  if(appdata$freqfilter == "Present"){
    snpdata = snpdata[snpdata$MajFreq > appdata$fixfreq,]  
  }
  if(appdata$freqfilter == "Start Polymorphic"){
    snpdata = snpdata[(snpdata$MajFreq - snpdata$MajDiff) > appdata$fixfreq & (snpdata$MajFreq - snpdata$MajDiff) < (1.0-appdata$fixfreq),]  
  }
  if(appdata$freqfilter == "Start Present"){
    snpdata = snpdata[(snpdata$MajFreq - snpdata$MajDiff) > appdata$fixfreq,]  
  }
  
  desc = c(appdata$freqfilter," ")

  ### Filter by direction
  if(appdata$plotpos == FALSE){ 
    snpdata = snpdata[snpdata$MajDiff <= 0,]
  }
  if(appdata$plotneg == FALSE){ 
    snpdata = snpdata[snpdata$MajDiff >= 0,]
  }
  if(appdata$plotfix == FALSE){ 
    snpdata = snpdata[snpdata$MajDiff != 0,]
  }
  desc =c(desc,paste0(c("|","+","-","=","|")[c(TRUE,appdata$plotpos,appdata$plotneg,appdata$plotfix,TRUE)],collapse=""))
  
  #1# Parents
  allparents = c("mbg344","mbg461","mbg474","mbg475","mbg479","mbg481","mbg482","mbg541","mbg542","mbg549","mbg557","mbg558","mbg602")
  desc = c(desc,"Par:",appdata$parlist)
  if(substr(appdata$parlist,1,3) == 'mbg'){   # Single Parent
    if(appdata$paruniq){  # Unique to parent
      if(appdata$parinv){
        snpdata = snpdata[snpdata$Parents != appdata$parlist,]  
        desc = c(desc," (UniqInv)")
      }else{
        snpdata = snpdata[snpdata$Parents == appdata$parlist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in parent
      if(appdata$parinv){
        snpdata = snpdata[regexpr(appdata$parlist,snpdata$Parents) < 0,]  
        desc = c(desc," (Inv)")
      }else{
        snpdata = snpdata[regexpr(appdata$parlist,snpdata$Parents) > 0,]  
      }
    }
  }
  
  if(appdata$parlist == "Any"){   # Any parent
    if(appdata$paruniq){  # Unique to parent
      if(appdata$parinv){
        snpdata = snpdata[! snpdata$Parents %in% allparents,]  
      }else{
        snpdata = snpdata[snpdata$Parents %in% allparents,]  
      }
    }
    # Present in 1+ parents
    else{
      if(appdata$parinv){
        snpdata = snpdata[regexpr("mbg",snpdata$Parents) < 0,]  
      }else{
        snpdata = snpdata[regexpr("mbg",snpdata$Parents) > 0,]  
      }
    }
  }  
  
  if(appdata$parlist == "All"){   # All parents
    for(mypar in allparents){
      if(appdata$parinv){
        snpdata = snpdata[regexpr(mypar,snpdata$Parents) < 0,]  
      }else{
        snpdata = snpdata[regexpr(mypar,snpdata$Parents) > 0,]  
      }
    }
  }  
  if(appdata$parlist == "None"){  # No parent
    snpdata = snpdata[snpdata$Parents == "",]  
  }
  
  
  #2# MBG non-parent Strains
  desc = c(desc,"; MBG:",appdata$strlist)
  if(substr(appdata$strlist,1,3) == 'mbg' & ! appdata$strcol){   # Single strain
    if(appdata$struniq){  # Unique to strain
      if(appdata$strinv){
        snpdata = snpdata[snpdata$Strains != appdata$strlist,]  
        desc = c(desc," (UniqInv)")
      }else{
        snpdata = snpdata[snpdata$Strains == appdata$strlist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in strain
      if(appdata$strinv){
        snpdata = snpdata[regexpr(appdata$strlist,snpdata$Strains) < 0,]  
        desc = c(desc," (Inv)")
      }else{
        snpdata = snpdata[regexpr(appdata$strlist,snpdata$Strains) > 0,]  
      }
    }
  }
  if(appdata$strlist == "Any" & ! appdata$strcol){   # Any strain
    if(appdata$strinv){
      # Set all rows to TRUE, then remove SNPs in strains
      keepers = regexpr("mbg",snpdata$Strains) > -2
      for(mypar in appdata$mbgstrainlist){
        keepers = keepers & regexpr(mypar,snpdata$Strains) < 0
      }
      snpdata = snpdata[keepers,]   
    }else{
      # Set all rows to FALSE, then add SNPs in strains
      keepers = regexpr("mbg",snpdata$Strains) < -2
      for(mypar in appdata$mbgstrainlist){
        keepers = keepers | regexpr(mypar,snpdata$Strains) > 0
      }
      snpdata = snpdata[keepers,]  
    }
  }  
  if(appdata$strlist == "All" & ! appdata$strcol){   # All strains
    writeLines(paste(appdata$mbgstrainlist))
    for(mypar in appdata$mbgstrainlist){
      if(appdata$strinv){
        snpdata = snpdata[regexpr(mypar,snpdata$Strains) < 0,]  
      }else{
        snpdata = snpdata[regexpr(mypar,snpdata$Strains) > 0,]  
      }
    }
  }
  if(appdata$strlist == "None" & ! appdata$strcol){  # No strain
    snpdata = snpdata[snpdata$Strains == "",]  
  }
  
  #2.2# MBG non-parent Strains 2
  desc = c(desc,"; MBG:",appdata$str2list)
  if(substr(appdata$str2list,1,3) == 'mbg'){   # Single strain
    if(appdata$str2uniq){  # Unique to strain
      if(appdata$str2inv){
        snpdata = snpdata[snpdata$Strains != appdata$str2list,]  
        desc = c(desc," (UniqInv)")
      }else{
        snpdata = snpdata[snpdata$Strains == appdata$str2list,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in strain
      if(appdata$str2inv){
        snpdata = snpdata[regexpr(appdata$str2list,snpdata$Strains) < 0,]  
        desc = c(desc," (Inv)")
      }else{
        snpdata = snpdata[regexpr(appdata$str2list,snpdata$Strains) > 0,]  
      }
    }
  }
  if(appdata$str2list == "Any"){   # Any strain
    if(appdata$str2inv){
      keepers = regexpr("mbg",snpdata$Strains) > -2
      for(mypar in appdata$mbgstrain2list){
        keepers = keepers & regexpr(mypar,snpdata$Strains) < 0
      }
      snpdata = snpdata[keepers,]   
    }else{
      keepers = regexpr("mbg",snpdata$Strains) < -2
      for(mypar in appdata$mbgstrain2list){
        keepers = keepers | regexpr(mypar,snpdata$Strains) > 0
      }
      snpdata = snpdata[keepers,]  
    }
  }  
  if(appdata$str2list == "All"){   # All strains
    writeLines(paste(appdata$mbgstrainlist))
    for(mypar in appdata$mbgstrainlist){
      if(appdata$str2inv){
        snpdata = snpdata[regexpr(mypar,snpdata$Strains) < 0,]  
      }else{
        snpdata = snpdata[regexpr(mypar,snpdata$Strains) > 0,]  
      }
    }
  }
  if(appdata$str2list == "None"){  # No strain
    snpdata = snpdata[snpdata$Strains == "",]  
  }

    #3# Populations
  desc = c(desc,"; Pop:",appdata$poplist)
  if(substr(appdata$poplist,1,3) == 'Pop'){   # Single population
    if(appdata$popuniq){  # Unique to population
      if(appdata$popinv){
        desc = c(desc," (UniqInv)")
        snpdata = snpdata[snpdata$Pops != appdata$poplist,]  
      }else{
        snpdata = snpdata[snpdata$Pops == appdata$poplist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in population
      if(appdata$popinv){
        snpdata = snpdata[regexpr(appdata$poplist,snpdata$Pops) < 0,]  
        desc = c(desc," (Inv)")
      }else{
        snpdata = snpdata[regexpr(appdata$poplist,snpdata$Pops) > 0,]  
      }
    }
  }
  if(appdata$poplist == "Any"){   # Any population
    if(appdata$popinv){
      snpdata = snpdata[regexpr("Pop",snpdata$Pops) < 0,]  
    }else{
      snpdata = snpdata[regexpr("Pop",snpdata$Pops) > 0,]  
    }
  }  
  if(appdata$poplist == "All"){   # All strains
    for(mypar in c("Pop00","Pop04","Pop06","Pop09","Pop10","Pop11","Pop12","Pop13")){
      if(appdata$popinv){
        snpdata = snpdata[regexpr(mypar,snpdata$Pops) < 0,]  
      }else{
        snpdata = snpdata[regexpr(mypar,snpdata$Pops) > 0,]  
      }
    }
  }
  if(appdata$poplist == "None"){  # No population
    if(appdata$popinv){
      snpdata = snpdata[regexpr("Pop",snpdata$Pops) > 0,]  
    }else{
      snpdata = snpdata[regexpr("Pop",snpdata$Pops) < 0,]  
    }
  }
  
  #4# SNP Types and Effects  
  if(appdata$snplist != "All" & ! appdata$plotcol %in% c("SNPType Selection")){
    if(appdata$snplist == "CDS"){
      snpdata = snpdata[snpdata$SNPType != "SNP",]  
    }else{
      snpdata = snpdata[snpdata$SNPType == appdata$snplist,]  
    }
  }
  if(appdata$efflist != "All" & ! appdata$plotcol %in% c("SNPEffect Selection")){
    snpdata = snpdata[snpdata$SNPEffect == appdata$efflist,]  
  }
  pop1 = colnames(snpdata)[6]
  pop2 = colnames(snpdata)[8]
  desc = paste0(c(desc,"; SNP:(",appdata$snplist,"|",appdata$efflist,")"),collapse="")
  
  chromdata = appdata$feattable
  xmin = appdata$xmin
  xmax = appdata$xmax
  ymin = appdata$ymin
  ymax = appdata$ymax
  
  # Locus_tag restriction
  #writeLines(appdata$findft)
  if(appdata$findft != ""){
    desc = paste(desc,"(Features)")
    writeLines("FilterSNPs locus_tag restriction")
    writeLines(appdata$findft)
    writeLines(">>>")
    if(regexpr(',',appdata$findft) > 0){
      # Multiple (new), e.g. "YHL049C,YHL026C,YOL051W"
      ftdata = snpdata[snpdata$Locus == "BLANK" & snpdata$Pos >= 1 & snpdata$Pos <= -1,]
      for(locustag in strsplit(appdata$findft,",",fixed=TRUE)[[1]]){
        writeLines(paste0("Looking for ",locustag))
        if(locustag %in% appdata$feattable$locus_tag){
          ftchrom = as.character(appdata$feattable[appdata$feattable$locus_tag == locustag,]$Chrom[1])
          fmin = min(appdata$feattable[appdata$feattable$locus_tag == locustag,]$Start)
          fmax = max(appdata$feattable[appdata$feattable$locus_tag == locustag,]$End)
          if(ftchrom %in% snpdata$Locus){
            ftdata = rbind(ftdata,snpdata[snpdata$Locus == ftchrom & snpdata$Pos >= fmin & snpdata$Pos <= fmax,])
          }
        }
      }
      if(nrow(ftdata) < 1){
        return(list(desc=paste("No selected features found in filtered SNP data."), snpdata=ftdata))
      }
      print(dim(ftdata))
      snpdata = ftdata          
      
    }else{
      # Single (old)
      writeLines("Single feature")
      ftchrom = as.character(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$Chrom[1])
      #writeLines(ftchrom)
      fmin = min(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$Start)
      #writeLines(fmin)
      fmax = max(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$End)
      writeLines("Looking for ftsnp")
      if(ftchrom %in% snpdata$Locus){
        snpdata = snpdata[snpdata$Locus == ftchrom & snpdata$Pos >= fmin & snpdata$Pos <= fmax,]  
        if(nrow(snpdata) < 1){
          return(list(desc=paste(ftchrom,"not found in filtered SNP data."), snpdata=snpdata))
        }
      }
      #else{
      #  return(paste(ftchrom,"not found in filtered SNP data."))
      #}
    }
  }
  
  desc = paste0(desc," ",pop1," to ",pop2," (n = ",nrow(snpdata),")")
  writeLines(desc)

  return( list(desc=desc, snpdata=snpdata) )
  #appdata$intro = desc
  #return(snpdata)
}



codedump = function(){
  #2# MBG non-parent Strains
  desc = c(desc,"; MBG:",appdata$strlist)
  if(substr(appdata$strlist,1,3) == 'mbg' & ! appdata$strcol){   # Single strain
    if(appdata$struniq){  # Unique to strain
      if(appdata$strinv){
        snpdata = snpdata[snpdata$Strains != appdata$strlist,]  
        desc = c(desc," (UniqInv)")
      }else{
        snpdata = snpdata[snpdata$Strains == appdata$strlist,]  
        desc = c(desc," (Uniq)")
      }
    }else{   # All SNPs in strain
      if(appdata$strinv){
        snpdata = snpdata[regexpr(appdata$strlist,snpdata$Strains) < 0,]  
        desc = c(desc," (Inv)")
      }else{
        snpdata = snpdata[regexpr(appdata$strlist,snpdata$Strains) > 0,]  
      }
    }
  }
  if(appdata$strlist == "Any" & ! appdata$strcol){   # Any strain
    if(appdata$strinv){
      # Set all rows to TRUE, then remove SNPs in strains
      keepers = regexpr("mbg",snpdata$Strains) > -2
      for(mypar in appdata$mbgstrainlist){
        keepers = keepers & regexpr(mypar,snpdata$Strains) < 0
      }
      if(appdata$struniq){  # Not uniquely present in any strains of interest
        # Cycle back through and add back SNPs in that strain AND OTHER strains too
        for(mypar in appdata$mbgstrainlist){
          others = appdata$mbgstrainlist[appdata$mbgstrainlist != mpar]
          for(notme in others){
            keepers = keepers | (regexpr(mypar,snpdata$Strains) > 0 & regexpr(notme,snpdata$Strains) > 0)
          }
        }
      }
      
      snpdata = snpdata[keepers,]   
    }else{
      # Set all rows to FALSE, then add SNPs in strains
      keepers = regexpr("mbg",snpdata$Strains) < -2
      for(mypar in appdata$mbgstrainlist){
        keepers = keepers | regexpr(mypar,snpdata$Strains) > 0
      }
      
      if(appdata$struniq){  # Uniquely present in any strains of interest
        # Cycle back through and remove SNPs in OTHER strains too
        for(mypar in appdata$mbgstrainlist){
          others = appdata$mbgstrainlist[appdata$mbgstrainlist != mpar]
          for(notme in others){
            keepers = keepers & regexpr(notme,snpdata$Strains) < 0
          }
        }
      }
      snpdata = snpdata[keepers,]  
    }
  }  
  if(appdata$strlist == "All" & ! appdata$strcol){   # All strains
    writeLines(paste(appdata$mbgstrainlist))
    for(mypar in appdata$mbgstrainlist){
      if(appdata$strinv){
        snpdata = snpdata[regexpr(mypar,snpdata$Strains) < 0,]  
      }else{
        snpdata = snpdata[regexpr(mypar,snpdata$Strains) > 0,]  
      }
    }
  }
  if(appdata$strlist == "None" & ! appdata$strcol){  # No strain
    snpdata = snpdata[snpdata$Strains == "",]  
  }
}


