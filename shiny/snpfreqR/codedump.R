# This is the filtering code removed from the plot code in place of snpfilter. Note that it built the description too, so might need reinstating. #


#>> FILTER START >>#
#i# This filtering is used for the data plots. Filtering for the SNP table will need to be updated separately #i#
### Subsetting data
if(appdata$freqfilter == "Fixed"){
  snpdata = snpdata[snpdata$MajFreq >= (1.0-appdata$fixfreq),]  
}
if(appdata$freqfilter == "Lost"){
  snpdata = snpdata[snpdata$MajFreq <= appdata$fixfreq,]  
}
if(appdata$freqfilter == "Polymorphic"){
  snpdata = snpdata[snpdata$MajFreq >= appdata$fixfreq & snpdata$MajFreq <= (1.0-appdata$fixfreq),]  
}
#1# Parents
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
  if(appdata$parinv){
    snpdata = snpdata[regexpr("mbg",snpdata$Parents) < 0,]  
  }else{
    snpdata = snpdata[regexpr("mbg",snpdata$Parents) > 0,]  
  }
}  
if(appdata$parlist == "None"){  # No parent
  snpdata = snpdata[snpdata$Parents == "",]  
}


#2# MBG non-parent Strains
desc = c(desc,"; MBG:",appdata$strlist)
if(substr(appdata$strlist,1,3) == 'mbg'){   # Single strain
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
      desc = c(desc," (nv)")
    }else{
      snpdata = snpdata[regexpr(appdata$strlist,snpdata$Strains) > 0,]  
    }
  }
}
if(appdata$strlist == "Any"){   # Any strain
  if(appdata$strinv){
    snpdata = snpdata[regexpr("mbg",snpdata$Strains) < 0,]  
  }else{
    snpdata = snpdata[regexpr("mbg",snpdata$Strains) > 0,]  
  }
}  
if(appdata$strlist == "None"){  # No strain
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
if(appdata$poplist == "None"){  # No population
  snpdata = snpdata[snpdata$Pops == "",]  
}

#4# SNP Types and Effects  
if(appdata$snplist != "All"){
  if(appdata$snplist == "CDS"){
    snpdata = snpdata[snpdata$SNPType != "SNP",]  
  }else{
    snpdata = snpdata[snpdata$SNPType == appdata$snplist,]  
  }
}
if(appdata$efflist != "All"){
  snpdata = snpdata[snpdata$SNPEffect == appdata$efflist,]  
}
pop1 = colnames(snpdata)[6]
pop2 = colnames(snpdata)[8]
desc = paste0(c(desc,"; SNP:(",appdata$snplist,"|",appdata$efflist,")"),collapse="")
desc = paste0(desc," ",pop1," to ",pop2," (n = ",nrow(snpdata[snpdata$Locus==chrom,])," / ",nrow(snpdata),")")
writeLines(desc)

chromdata = appdata$feattable
comparison = paste0(appdata$complabel," ",pop1," to ",pop2," (n = ",nrow(snpdata[snpdata$Locus==chrom,])," / ",nrow(snpdata),")")
xmin = appdata$xmin
xmax = appdata$xmax
ymin = appdata$ymin
ymax = appdata$ymax

# Locus_tag restriction
if(appdata$findft != ""){
  writeLines("locus_tag restriction")
  writeLines(appdata$findft)
  writeLines(">>>")
  if(regexpr(',',appdata$findft) > 0){
    # Multiple (new), e.g. "YHL049C,YHL026C,YOL051W"
    ftdata = snpdata[snpdata$Locus == "BLANK" & snpdata$Pos >= 1 & snpdata$Pos <= -1,]
    for(locustag in strsplit(appdata$findft,",",fixed=TRUE)){
      writeLines(paste0("Looking for",locustag))
      ftchrom = as.character(appdata$feattable[appdata$feattable$locus_tag == locustag,]$Chrom[1])
      fmin = min(appdata$feattable[appdata$feattable$locus_tag == locustag,]$start)
      fmax = max(appdata$feattable[appdata$feattable$locus_tag == locustag,]$end)
      if(ftchrom %in% snpdata$Locus){
        ftdata = rbind(ftdata,snpdata[snpdata$Locus == ftchrom & snpdata$Pos >= fmin & snpdata$Pos <= fmax,])
        if(nrow(ftdata) < 1){
          return(paste(ftchrom,"not found in filtered SNP data."))
        }
      }else{
        return(paste(ftchrom,"not found in loaded data."))
      }
    }
    snpdata = ftdata          
    
  }else{
    # Single (old)
    ftchrom = as.character(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$Chrom[1])
    #writeLines(ftchrom)
    fmin = min(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$Start)
    #writeLines(fmin)
    fmax = max(appdata$feattable[appdata$feattable$locus_tag == appdata$findft,]$End)
    writeLines("Looking for ftsnp")
    if(ftchrom %in% snpdata$Locus){
      snpdata = snpdata[snpdata$Locus == ftchrom & snpdata$Pos >= fmin & snpdata$Pos <= fmax,]  
      if(nrow(snpdata) < 1){
        return(paste(ftchrom,"not found in filtered SNP data."))
      }
    }else{
      return(paste(ftchrom,"not found in filtered SNP data."))
    }
  }
}
#<< FILTER END <<#