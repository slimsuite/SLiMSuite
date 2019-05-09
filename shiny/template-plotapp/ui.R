#i# See the `main.R` file for version, contact and license information.
#i# The `main.R` script loads libraries and contains the initial parameter settings and functions.
source("main.R")

############### ::: USER INTERFACE ::: ##################
# Define the UI.
shinyUI(fluidPage(
  #>>># This code block should be copied to the standalone app ui code when changed #>>>#
  # Application title.
  titlePanel(info$apptitle),
  # Sidebar layout with parameter input objects and tabs of output
  sidebarLayout(

    # Sidebar with parameter input objects
    sidebarPanel(
      # Description markdown file
      #i# Can include static MD like this: #includeMarkdown("summary.md"),
      htmlOutput("summary"),
      selectInput("freqlist", "Timepoint:", "Load data!", width='300px'),#multiple=TRUE),
      selectInput("chromlist", "Chromosome:", c(settings$chrom), width='300px'),
      
      #radioButtons("loadplot", "Load or Plot data", choices = c("Load","Plot"), selected = "Load", inline = TRUE, width = NULL, choiceNames = NULL, choiceValues = NULL),

      # These parameters alter plot appearance but do not change the data
      checkboxInput("plotdata", "Show plot/filters (toggle with Load options)", value=FALSE),
      
      #i# Load Data
      conditionalPanel(
        condition = "input.plotdata == false",
        checkboxInput("multiload", "Load Multiple SNPFreq tables", value=TRUE),

        conditionalPanel(
          condition = "input.multiload == false",
          
          fileInput("snpfreq","SNP Frequency Table"),#multiple=TRUE),
          fileInput("features","Features Table")
          #x# No need for this?:
          #x# fileInput("snpmap","Full SNPMap Table")
        ),
        
        actionButton("load", "Load Data")
      ),
      
      conditionalPanel(
        condition = "input.plotdata == true",
        #actionButton("makeplot", "Plot data (Plot tab)"),
        
        
        # Plot settings
        htmlOutput("hr"),
        checkboxInput("plotopt", "Show plot options", value=TRUE),
        checkboxInput("plotzoom", "Show zoom settings", value=FALSE),
        checkboxInput("plotfilt", "Show SNP filtering settings", value=FALSE),
        conditionalPanel(
          condition = "input.plotopt == true",
          htmlOutput("plotoptmd"),
          textInput("title", "Plot title:", "SNPFreq"),
          selectInput("titlestyle", "Plot title style:", c("Standard","Pure","Info","None"), "Info", width='300px'),
          #selectInput("plottype", "Plot strategy:", c("Standard","Time Lapse","Chrom Cycle","All Chrom (PNG Only)","Monster (PNG Only)","Genome multi-plot","Time multi-plot"), width='300px'),   #i# ,"Time Lapse" does not work!
          selectInput("plottype", "Plot strategy:", c("Standard","Time Lapse","Chrom Cycle","All Chrom (PNG Only)","Monster (PNG Only)"), width='300px'),   #i# ,"Time Lapse" does not work!
          selectInput("plotfield", "Plot data:", c("SNPFreq","Starting SNPFreq","Ending SNPFreq","SNPFreq Change","Absolute SNPFreq Change"), width='300px'),   #i# ,"Time Lapse" does not work!
          selectInput("plothist", "Plot type:", c("Scatterplot","Histogram","Histogram (all chrom)","TimeLine","TimeLine (all chrom)"), width='300px'),   #i# ,"Time Lapse" does not work!
          selectInput("plotcol", "Plotting colours:", c("Direction","Parent","MBG Strain","Fixation","Strain 1 Selection","No Parent","SNPType Selection","SNPEffect Selection"), "Direction",width='300px'),
          checkboxInput("plotpos", "Plot positive freq changes", value=settings$plotpos),
          checkboxInput("plotneg", "Plot negative freq changes", value=settings$plotneg),
          checkboxInput("plotfix", "Plot unchanged freq positions", value=settings$plotneg),
          checkboxInput("plothighlights", "Backfill allele markers", value=TRUE),
          checkboxInput("pauseplot", "Pause plotting (whilst changing options)", value=FALSE)
        ),
        conditionalPanel(
          condition = "input.plotzoom == true",
          htmlOutput("zoomoptmd"),
          numericInput("xmin", "Chromosome start pos (kb):", settings$xmin, min = 0, max = 1000, step = 1),
          numericInput("xmax", "Chromosome end pos (kb) [-1 for end]:", settings$xmax, min = -1, max = 100, step = 1)
        ),
        conditionalPanel(
          condition = "input.plotfilt == true",
          htmlOutput("filtoptmd"),
          selectInput("poplist", "Restrict to Populations:", c("No Filter","All","Any","None","Pop00","Pop04","Pop06","Pop09","Pop10","Pop11","Pop12","Pop13"), "Any",width='200px'),
          checkboxInput("popuniq", "Restrict to population-unique SNPs", value=FALSE),
          checkboxInput("popinv", "Invert Population selection", value=FALSE),
          selectInput("parlist", "Restrict to Parents:", c("No Filter","All","Any","None","mbg344","mbg461","mbg474","mbg475","mbg479","mbg481","mbg482","mbg541","mbg542","mbg549","mbg557","mbg558","mbg602"), "No Filter",width='200px'),
          checkboxInput("paruniq", "Restrict to parent-unique SNPs", value=FALSE),
          checkboxInput("parinv", "Invert Parent selection", value=FALSE),
          #selectInput("strlist", "Restrict to MBG Strains:", c("All","Any","None",'mbg11a','mbg1871','mbg2303','mbgag26','mbgag35','mbgh207'), "All",width='200px'),
          selectInput("strlist", "Restrict to MBG Strains (1):", c("No Filter","All","Any",strsplit(settings$mbgstrains,',',TRUE)[[1]]), "No Filter",width='200px'),
          checkboxInput("struniq", "Restrict to strain-unique SNPs (not limited to list)", value=FALSE),
          checkboxInput("strinv", "Invert Strain selection", value=FALSE),
          checkboxInput("strcol", "Colour Strain selection only (single strains)", value=FALSE),
          selectInput("str2list", "Restrict to MBG Strains (2):", c("No Filter","All","Any",'mbg11a','mbg1871','mbg2303'), "No Filter",width='200px'),
          checkboxInput("str2uniq", "Restrict to strain-unique SNPs (not limited to list)", value=FALSE),
          checkboxInput("str2inv", "Invert Strain selection", value=FALSE),
          selectInput("snplist", "Restrict to SNPType:", c("All"), "All",width='200px'),
          selectInput("efflist", "Restrict to SNPEffect:", c("All","gene","mRNA","CDS","rRNA","tRNA","ncRNA","mobile","LTR","origin","centromere","telomere"), "No Filter",width='200px'),
          #numericInput("fixfreq", "Frequency Buffer for fixation",0.05,min=0,max=1,step=0.005),
          #Moved to plot tab: selectInput("freqfilter", "Final freq for SNPs to plot:", c("All","Fixed","Lost","Polymorphic"),"All",width='200px'),
          textInput("findft", "Restrict to Feature (locus_tag):", "")
        )
        
      ),

      
      # HTML content to include at end. Contains version number.
      htmlOutput("footer")
      #i# Can include static HTML like this: includeHTML("footer.html")
    ),
    
    # Main output panel
    mainPanel(
      # Overall run info (data$intro)
      verbatimTextOutput("intro"),
      # Output tabs
      tabsetPanel(
        # Markdown descriptions of tool
        tabPanel("Instructions", 
                 checkboxInput("showdesc", "Show App instructions", value=TRUE),
                 checkboxInput("showinfo", "Show explanation of App output", value=FALSE),
                 conditionalPanel(
                   condition = "input.showdesc == true",
                   includeMarkdown("description.md"),
                   htmlOutput("debug")
                 ),
                 
                 conditionalPanel(
                   condition = "input.showinfo == true",
                   includeMarkdown("information.md")      
                   #htmlOutput("outfmt")
                 )
        ),
        # Main Plot Panel
        tabPanel("Plot", 
           conditionalPanel(
             condition = "input.plotfilt == true",
             selectInput("freqfilter", "Final freq for SNPs to plot:", c("All","Fixed","Lost","Polymorphic","Present","Start Polymorphic","Start Present"),"All",width='300px'),
             numericInput("fixfreq", "Frequency Buffer for fixation",0.05,min=0,max=1,step=0.005)
           ),
           conditionalPanel(
             condition = "input.plotdata == true",
             actionButton("makeplot", "Plot data (Plot tab)"),
             htmlOutput("hr5"),
             plotOutput("freqplot",height=settings$plotheight),
             htmlOutput("hr6"),
             actionButton("savepng", "Save plot as PNG (See Adv. Settings)")
          )
        ),
        tabPanel("SNP Frequencies", 
                 checkboxInput("filterfreqtable", "Filter SNPs in table", value=TRUE),
                 checkboxInput("zoomfreqtable", "Restrict table to chromosome and zoom region", value=TRUE),
                 dataTableOutput("freqtable")
        ),
        tabPanel("Timeline", 
           checkboxInput("filtertimetable", "Filter SNPs in table", value=TRUE),
           checkboxInput("zoomtimetable", "Restrict table to chromosome and zoom region", value=TRUE),
           dataTableOutput("timetable")
        ),
        tabPanel("Features", dataTableOutput("feattable")),
        tabPanel("Log", verbatimTextOutput("activitylog")),
        #tabPanel("SNPMap", dataTableOutput("snptable")),
        tabPanel("Adv. Settings", 
                 # Input path settings
                 htmlOutput("pathset"),
                 conditionalPanel(
                   condition = "input.multiload == true",
                   #i# Edit this during testing. Current default path is set for Google Drive sharing
                   textInput("freqpath", "Path to SNPFreq files:", settings$inputpath, width='600px'),
                   textInput("freqbase", "Root filename for SNPFreq files:", "pop00Freq", width='300px')
                   #x#shinyDirButton("freqpath", "Chose directory", "Upload")
                 ),
                 textInput("mbgstrainlist", "MBG strains (1):", settings$mbgstrains, width='600px'),
                 textInput("mbgstrain2list", "MBG strains (2):", settings$mbgstrains, width='600px'),
                 #textInput("multipoplist", "Population comparisons for multipanel population plots:", paste0(settings$popcomb,collapse=","), width='600px'),
                 #textInput("multichromlist", "Chromosome prefixes for multipanel chromosome plots:", paste0(settings$multichrom,collapse=","), width='600px'),
                 
                 # PNG settings
                 htmlOutput("hr2"),
                 checkboxInput("plotpng", "Show PNG output settings", value=TRUE, width='600px'),
                 
                 conditionalPanel(
                   condition = "input.plotpng == true",
                   textInput("pngpath", "PNG Path:", settings$pngpath, width='600px'),
                   textInput("pngbase", "PNG Basefile:", "SNPFreq", width='600px'),
                   numericInput("pngheight", "Height (px) for PNG:", settings$pngheight, min = 400, max = 2400, step = 100),
                   numericInput("pngwidth", "Width (px) for PNG:", settings$pngwidth, min = 400, max = 3600, step = 100),
                   numericInput("pointsize", "Point size for PNG text:", settings$pointsize, min = 1, max = 73, step = 1)
                 ),
                 
                 # Plot settings
                 htmlOutput("hr3"),
                 checkboxInput("plotadv", "Show adv. plotting settings", value=FALSE),
                 conditionalPanel(
                   condition = "input.plotadv == true",
                    #numericInput("plotheight", "Plot height (px):", settings$plotheight, min = 100, max = 4800, step = 100),
                    numericInput("pcut", "Min change probability for shading:", settings$pcut, min = 0, max = 1, step = 0.01),
                    numericInput("ymin", "Min Freq (y-axis):", settings$ymin, min = 0, max = 1, step = 0.1),
                    numericInput("ymax", "Max Freq (y-axis):", settings$ymax, min = 0, max = 1, step = 0.1)
                 )
        )
      )
    )
  )
  #<<<# End of ui.R code block #<<<#
))


############### ::: SERVER CODE ::: ##################
#i# server.R defines server logic

