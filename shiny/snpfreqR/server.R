#i# See the `main.R` file for version, contact and license information.
#i# The `main.R` script loads libraries and contains the initial parameter settings and functions.
options(shiny.maxRequestSize=50*1024^2)
source("main.R")

# This is the code that goes inside the server object
shinyServer(function(input, output, session) {

  ### SECTION 1 - Setup a reactive data list for the server data
  adata <- reactiveValues(
      data = setupData()
  )
  #i# Can update in a render function using adata$data = ...

  # HEADER: Verbatim text outputs
  output$intro <- renderText({

    adata$data$intro
    input$load
    if(input$load > adata$data$loads & input$plotdata == FALSE){
      isolate({
        withProgress(message="Loading data", value=0, {
          if(input$multiload == FALSE){
            writeLines(paste("Loading",nrow(input$snpfreq) ,"data",input$load))
            #i# Load SNP Frequencies
            freqdb = read.delim(input$snpfreq$datapath,header=TRUE,stringsAsFactors = TRUE,sep="\t")
            # if(! 'AltPos' %in% colnames(freqdb)){
            #   freqdb$AltPos = freqdb$Pos
            #   freqdb$AltLocus = freqdb$Locus
            # }
            # freqdb$Pops = as.character(freqdb$Pops)
            # if(length(freqdb[freqdb$Pops == "Pops",]$Pops) > 0){
            #   freqdb[freqdb$Pops == "Pops",]$Pops = ""
            # }
            # freqdb$Pops = as.factor(freqdb$Pops)
            # freqdb$SNPEffect = as.character(freqdb$SNPEffect)
            # if(length(freqdb[freqdb$SNPType == "NS",]$SNPEffect) > 0){
            #   freqdb[freqdb$SNPType == "NS",]$SNPEffect = "Nonsyn"
            # }
            # if(length(freqdb[freqdb$SNPType == "SYN",]$SNPEffect) > 0){
            #   freqdb[freqdb$SNPType == "SYN",]$SNPEffect = "Syn"
            # }
            # if(length(freqdb[freqdb$SNPType == "NON",]$SNPEffect) > 0){
            #   freqdb[freqdb$SNPType == "NON",]$SNPEffect = "STOP"
            # }
            # if(length(freqdb[freqdb$SNPType == "EXT",]$SNPEffect) > 0){
            #   freqdb[freqdb$SNPType == "EXT",]$SNPEffect = "Extension"
            # }
            # freqdb$SNPEffect = as.factor(freqdb$SNPEffect)
            freqdb = freqTableCleanup(freqdb)
            adata$data$freqloaded = TRUE
            adata$data$freqtable = freqdb
            writeLines(input$snpfreq$name)
            incProgress(1/3)
            #i# Load Features Table
            ftdb = read.delim(input$features$datapath,header=TRUE,stringsAsFactors = TRUE,sep="\t")
            writeLines(input$features$name)
            updateSelectInput(session, "freqlist",
              choices = c("Loaded data"),
              selected = "Loaded data"
            )            
            adata$data$freqdb = list("Loaded data"=freqdb)
            incProgress(1/3)
          }else{
            #x#shinyDirChoose(input, 'freqpath', roots = c(home = '~'), filetypes = c('', 'txt', 'tdt'))
            writeLines(paste("Loading",input$freqpath,"data",input$load))
            ftfile = paste(input$freqpath,"ref.Feature.tdt",sep="/")
            ftfile = paste(input$freqpath,"sgd.R64.2.1.Feature.tdt",sep="/")
            ftdb = read.delim(ftfile,header=TRUE,stringsAsFactors = TRUE,sep="\t")
            writeLines(ftfile)
            incProgress(1/3)
            freqfile = paste(input$freqpath,"allfreq.Pop00.vs.Pop04.snpfreq.tdt",sep="/")
            freqdb = read.delim(freqfile,header=TRUE,stringsAsFactors = TRUE,sep="\t")
            freqdb = freqTableCleanup(freqdb)
            adata$data$freqloaded = TRUE
            adata$data$freqtable = freqdb
            writeLines(freqfile)
            incProgress(1/21)
            
            popcomb = c("Pop00.vs.Pop04","Pop04.vs.Pop06","Pop06.vs.Pop09","Pop09.vs.Pop10","Pop10.vs.Pop11","Pop11.vs.Pop12","Pop12.vs.Pop13")
            adata$data$freqdb = list()
            adata$data$freqdb[["Pop00.vs.Pop04"]] = freqdb
            for(combo in popcomb[2:7]){
              freqfile = paste0(input$freqpath,"/allfreq.",combo,".snpfreq.tdt")
              freqdb = read.delim(freqfile,header=TRUE,stringsAsFactors = TRUE,sep="\t")
              freqdb = freqTableCleanup(freqdb)
              adata$data$freqdb[[combo]] = freqdb  
              writeLines(freqfile)
              incProgress(1/21)
            }
            # freqname = as.character(input$snpfreq$name)
            # adata$data$freqfiles = c(adata$data$freqfiles,freqname)
            # print(adata$data$freqfiles)
            # adata$data$freqdb[[freqname]] = freqdb
            updateSelectInput(session, "freqlist",
                              #                  label = "Chromosome to Plot:",
                              choices = popcomb,
                              selected = "Pop00.vs.Pop04"
            )
          }

          #selectInput("snplist", "Restrict to SNPType:", c("All"), "text"),
          updateSelectInput(session, "snplist",
          #                  label = "Chromosome to Plot:",
                              choices = c("All","CDS",levels(freqdb$SNPType)),
                              selected = "All"
          )
          #selectInput("efflist", "Restrict to SNPEffect:", c("All"), "text"),
          updateSelectInput(session, "efflist",
          #                  label = "Chromosome to Plot:",
                              choices = c("All",levels(freqdb$SNPEffect)),
                              selected = "All"
          )
          
          # Update feature data                    
          colnames(ftdb) = c("Locus","feature","position","Start","End","locus_tag","protein_id","details")
          loci = as.character(levels(freqdb$Locus))
          accs = c()
          for(locus in loci){
            accs = c(accs,strsplit(locus,'_')[[1]][4])
          }
          ftdb$Chrom = as.character(ftdb$Locus)
          for(locus in as.character(levels(ftdb$Locus))){
            if(locus %in% accs){
              ftdb[ftdb$Locus == locus,]$Chrom = loci[which(accs==locus)]
            }
          }
          ftdb$Chrom = as.factor(ftdb$Chrom)
          ftdb$Source = "Ref"
          writeLines(as.character(dim(ftdb)))
          adata$data$altchr = levels(as.factor(as.character(ftdb[ftdb$Source=="Ref",]$Chrom)))
          adata$data$refchr = levels(as.factor(as.character(ftdb[ftdb$Source=="Ref",]$Chrom)))
          adata$data$featloaded = TRUE
          adata$data$feattable = ftdb
          incProgress(1/6)
          
          updateSelectInput(session, "chromlist",
                            label = "Chromosome to Plot:",
                            choices = c(levels(as.factor(as.character(ftdb[ftdb$Source=="Ref",]$Chrom))),"None"),
                            selected = settings$chrom
          )
          
          #i# Load Features Table
          if(is.null(input$snpmap) == FALSE){
            snpdb = read.delim(input$snpmap$datapath,header=TRUE,stringsAsFactors = TRUE,sep="\t")
            adata$data$snploaded = TRUE
            adata$data$snptable = snpdb
          }
          incProgress(1/6)
          adata$data$loads = input$load

          adata$data$intro = "Data Loaded. Ready to Plot Data."
        })  
      })
    }

    input$savepng
    if(input$savepng > adata$data$pngsaves & input$plotdata == TRUE){
      isolate({
        withProgress(message="Saving plot to PNG", value=0, {
          writeLines(paste(c("Saving plot to PNG",input$savepng)))
          adata$data$pngpath = input$pngpath
          adata$data$pngbase = input$pngbase
          adata$data$pngwidth = input$pngwidth
          adata$data$pngheight = input$pngheight
          adata$data$pointsize = input$pointsize
          if(input$plottype == "Time Lapse"){
            popcomb = c("Pop00.vs.Pop04","Pop04.vs.Pop06","Pop06.vs.Pop09","Pop09.vs.Pop10","Pop10.vs.Pop11","Pop11.vs.Pop12","Pop12.vs.Pop13")
            for(combo in popcomb){
              writeLines(combo)
              adata$data$freqlist = combo
              adata$data$intro = zoomPlot(adata$data,makepng = TRUE)
              incProgress(1/7)
            }
          }else{
            adata$data$intro = zoomPlot(adata$data,makepng = TRUE)
            incProgress(1/1)
          }
          adata$data$pngsaves = input$savepng
        })
      })
    }
          
    
    return(adata$data$intro)
  })
  
  
  
  
  output$freqtable = renderDataTable({
    # writeLines(adata$data$freqloaded)
    # print(input$freqlist)
      print(input$freqlist)
      adata$data$freqtable = adata$data$freqdb[[input$freqlist]]
      if(input$filterfreqtable){
        #Update plot settings
        adata$data$chrom = input$chromlist
        adata$data$plotpos = input$plotpos
        adata$data$plotneg = input$plotneg
        adata$data$plotfix = input$plotfix

        adata$data$freqfilter = input$freqfilter
        adata$data$fixfreq = input$fixfreq
        adata$data$parlist = input$parlist
        adata$data$strlist = input$strlist
        adata$data$poplist = input$poplist
        adata$data$paruniq = input$paruniq
        adata$data$struniq = input$struniq
        adata$data$popuniq = input$popuniq
        adata$data$parinv = input$parinv
        adata$data$strinv = input$strinv
        adata$data$popinv = input$popinv
        adata$data$snplist = input$snplist
        adata$data$efflist = input$efflist
        
        adata$data$findft = input$findft
        fdata = filterSNPs(adata$data,adata$data$freqtable)
      }else{
        fdata = adata$data$freqtable
      }
      if(input$zoomfreqtable){
        minpos = input$xmin
        maxpos = input$xmax
        if(maxpos < 1){ maxpos = max(fdata$Pos) }
        fdata = fdata[fdata$Pos >= minpos & fdata$Pos <= maxpos & fdata$Locus == input$chromlist,]
      }
      fdata
    },
    rownames=FALSE,
    options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 25)
  )
  
  output$feattable = renderDataTable({
    adata$data$feattable
  },
  rownames=FALSE,
  options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 25)
  )
  
  output$snptable = renderDataTable({
    adata$data$snptable
  },
  rownames=FALSE,
  options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 25)
  )

  PlotHeight = reactive(
    return( input$plotheight )
  )
  
  # Plot
  output$freqplot = renderPlot({
    #plot(c(0,1),c(0,1))
    input$makeplot
    if(input$makeplot > 0 & input$plotdata == TRUE){
      writeLines(input$chromlist)
      writeLines(input$titlestyle)
      input$freqfilter
      input$freqlist
      isolate({
  
        #!#adata$freqtable = adata$freqdb[[input$freqlist]]
        #Find feature
        adata$data$findft = ""
        if(input$findft != "" & input$findft %in% adata$data$feattable$locus_tag){
          ftchrom = as.character(adata$data$feattable[adata$data$feattable$locus_tag==input$findft,]$Chrom[1])
          writeLines(paste(input$findft,"=",ftchrom))
          if(input$chromlist != ftchrom){
            updateSelectInput(session, "chromlist",
                              selected = ftchrom
            )
          }
          adata$data$findft = input$findft
        }
        
        #Update plot settings
        adata$data$chrom = input$chromlist
        adata$data$titlestyle = input$titlestyle
        adata$data$xmin = input$xmin
        adata$data$xmax = input$xmax
        adata$data$plotpos = input$plotpos
        adata$data$plotneg = input$plotneg
        adata$data$plotfix = input$plotfix
        adata$data$plotheight = input$plotheight
        adata$data$ymin = input$ymin
        adata$data$ymax = input$ymax
        adata$data$complabel = input$title

        adata$data$freqfilter = input$freqfilter
        adata$data$fixfreq = input$fixfreq
        adata$data$parlist = input$parlist
        adata$data$strlist = input$strlist
        adata$data$poplist = input$poplist
        adata$data$paruniq = input$paruniq
        adata$data$struniq = input$struniq
        adata$data$popuniq = input$popuniq
        adata$data$parinv = input$parinv
        adata$data$strinv = input$strinv
        adata$data$popinv = input$popinv
        adata$data$snplist = input$snplist
        adata$data$efflist = input$efflist
        
        #Generate plot
        if(input$plottype == "Time Lapse" & input$freqlist != "Loaded data"){
          popcomb = c("Pop00.vs.Pop04","Pop04.vs.Pop06","Pop06.vs.Pop09","Pop09.vs.Pop10","Pop10.vs.Pop11","Pop11.vs.Pop12","Pop12.vs.Pop13")
          findcomb = which(popcomb == adata$data$freqlist)
          if(length(findcomb) < 1){
            combo = popcomb[1]
          }else{
            if(findcomb == 7){
              combo = popcomb[1]
            }else{
              combo = popcomb[findcomb+1]
            }
          }
          writeLines(combo)
          writeLines(adata$data$freqlist)
          adata$data$freqlist = combo
          adata$data$intro = zoomPlot(adata$data)
        }else{
          adata$data$freqlist = input$freqlist
          adata$data$intro = zoomPlot(adata$data)
          adata$data$multiplotx = 0
        }
        writeLines("Plot Rendered.")
        writeLines("")
      })
    }
  })
  
  # output$plot.ui <- renderUI({
  #   plotOutput("freqplot", height = PlotHeight())
  # })
  
  
  #i# Check whether a jobID looks legit and return True or False
  # isJobID <- function(jobid){
  #i# Check whether Job has run
  # checkJob <- function(jobid,password=""){
  #i# Function for returning the REST keys
  # getRestKeys <- function(jobid,password=""){
  #i# Return an R object with REST output
  # getRestOutput <- function(jobid,rest,outfmt="text",password=""){
    
  ### SECTION 2 - Status panel: response to Retrieve Button
  output$status <- renderText({
    input$retrieve
    if(input$retrieve > 0){
      isolate({
        withProgress(message="Checking JobID", value=0, {
          adata$data <- setupData()
          incProgress(1/4)
          #i# First, check whether it looks like a JobID
          if(isJobID(input$jobid) == FALSE){
            adata$data$status = paste("ERROR:",input$jobid,"is an invalid JobID.")
            return(paste(as.character(adata$data$status),sep="\n",collapse="\n"))
          }
          incProgress(1/4)
          #i# Next, check Job for completion
          jcheck = checkJob(input$jobid,input$password)
          if(jcheck != TRUE){
            adata$data$status = jcheck
            return(paste(as.character(adata$data$status),sep="\n",collapse="\n"))
          }
          incProgress(1/4)
          adata$data$restkeys = c(getRestKeys(input$jobid,input$password),settings$restkeys)
          incProgress(1/4)
        })  
        progx = length(adata$data$restkeys)
        withProgress(message="Retrieving data", value=0, {
          for(ikey in adata$data$restkeys){
            adata$data[[ikey]] = getRestOutput(input$jobid,ikey,password=input$password)
            incProgress(1/progx)
          }
        })
      })
    }
    return(paste(as.character(adata$data$status),sep="\n",collapse="\n"))
  })
  
  ### SECTION 3 - Output tabs: data rendering
  ### Standard Server Outputs
  output$ini <- renderText({
    inihead = c("ini file:","---------")
    return(paste(c(inihead,as.character(adata$data$ini)),sep="\n",collapse="\n")) 
  })
  
  output$outfmt <- renderUI({
    pretext = paste(as.character(adata$data$outfmt),sep="\n",collapse="\n")
    myhtml = c(paste0("<h2>",info$apptitle," Outputs</h2>"),
               paste0("<pre>",pretext,"</pre>")
    )
    return(HTML(paste(myhtml,sep="<br/>\n",collapse="<br/>\n")))
  })
  
  output$restkeys <- renderUI({
    pretext = paste(as.character(adata$data$restkeys),sep="\n",collapse="\n")
    myhtml = c(paste0("<h2>",adata$data$prog," Outputs</h2>"),   
               paste0("<pre>",pretext,"</pre>")
    )
    return(HTML(paste(myhtml,sep="<br/>\n",collapse="<br/>\n")))
  })
  
  output$help <- renderUI({
    pretext = paste(as.character(adata$data$help),sep="\n",collapse="\n")
    myhtml = c(paste0("<h2>",info$apptitle," Help</h2>"),
               paste0("<pre>",pretext,"</pre>")
    )
    return(HTML(paste(myhtml,sep="<br/>\n",collapse="<br/>\n")))
  })
  
  #i# Return multi-line HTML coded text
  output$wisdoms <- renderUI({
    myhtml = c("<h2>Random wisdom generator:</h2>\n<p>")
    myhtml = c(myhtml,paste(as.character(adata$data$wisdoms),sep="</p>\n<p>",collapse="</p>\n<p>","</p>\n"))
    HTML(paste0(myhtml))
    #return(HTML(paste(as.character(adata$data$wisdoms[,1]),sep="<br/>\n",collapse="<br/>\n")))
  })
  
  output$restout <- renderUI({
    myhtml = c(paste0("<h3>",input$restout," (",input$restformat,"):</h3>\n<p>"))
    myhtml = c(myhtml,paste(as.character(adata$data[[input$restout]]),sep="</p>\n<p>",collapse="</p>\n<p>","</p>\n"))
    HTML(paste0(myhtml))
    #return(HTML(paste(as.character(adata$data$wisdoms[,1]),sep="<br/>\n",collapse="<br/>\n")))
  })
  
  output$restoutTable = renderDataTable({
    adata$data[[input$restout]]
    },
    rownames=FALSE,
    options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 25)
  )
  
  output$warnings <- renderText({
    return(paste(as.character(adata$data$warnings),sep="\n",collapse="\n")) 
  })
  output$errors <- renderText({
    return(paste(as.character(adata$data$errors),sep="\n",collapse="\n")) 
  })
  
  # Text outputs
  output$errortext <- renderText({
    return("Error messages (if any):")
  })
  
  # Log output
  
  # If we ever need to make the data sortable etc.
  #i# tabPanel("Log", dataTableOutput("log")),
  output$log = renderDataTable({
    adata$data$log
    },
    rownames=FALSE,
    options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 25)
  )
  
  ### SECTION 4 - Footers and extra info
  output$footer <- renderUI({
    myhtml = c(paste("<hr>\n<p>&copy; 2017", info$author), # Richard J. Edwards",
               "An <a href=\"http://www.slimsuite.unsw.edu.au/shiny.php\">EdwardsLab</a> Shiny App",
               paste("Version", info$version, "| <a href=\"http://www.gnu.org/licenses/#GPL\">GNU GPL v3</a></p>"))
    return(HTML(paste(myhtml,sep="<br/>\n",collapse="<br/>\n")))
  })
  
  output$summary <- renderUI({
    return(HTML(renderMarkdown(text=info$description)))
  })
  
  output$plotoptmd <- renderUI({
    md = c("---","","#### Basic plot options:","")
    md = paste(md,sep="\n")
    return(HTML(renderMarkdown(text=md)))
  })

  
  output$zoomoptmd <- renderUI({
    md = c("---","","#### Plot/table zoom options:","")
    md = paste(md,sep="\n")
    return(HTML(renderMarkdown(text=md)))
  })

  output$filtoptmd <- renderUI({
    md = c("---","","#### SNP filtering options:","")
    md = paste(md,sep="\n")
    return(HTML(renderMarkdown(text=md)))
  })
  
  output$hr <- renderUI({
    return(HTML("<hr>"))
  })
  output$hr2 <- renderUI({ return(HTML("<hr>")) })
  output$hr3 <- renderUI({ return(HTML("<hr>")) })

  output$hr5 <- renderUI({ return(HTML("<hr>")) })
  output$hr6 <- renderUI({ return(HTML("<hr>")) })

    #i# Markdown description of setting the paths to the files
  output$pathset <- renderUI({
    pathmd = c("### Advanced settings","","The multiload function of SNPFreqR expects the `sgd.R64.2.1.Feature.tdt` and `allfreq.*.snpfreq.tdt` files to be present in a single directory. Set the path to these files below. Alternatively, uncheck the **Load Multiple SNPFreq tables** box to select a single timepoint and feature table.")
    pathmd = paste(pathmd,sep="\n")
    return(HTML(renderMarkdown(text=pathmd)))
  })
  
  
  #i# Debugging output
  output$debug <- renderUI({
    if(settings$debug == FALSE){ return("") }
    itxt = c()
    for(setkey in c("showdesc","showini","showout","showinfo")){
      itxt = c(itxt,paste0(setkey,"=",input[[setkey]]))
    }
    pretext = paste(itxt,sep="\n",collapse="\n")
    myhtml = c(paste0("<h3>Debugging</h3>"),
               paste0("<pre>",pretext,"</pre>")
    )
    return(HTML(paste(myhtml,sep="<br/>\n",collapse="<br/>\n")))
  })
  
    
  
})
