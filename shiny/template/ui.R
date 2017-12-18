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
      # Set parameters for REST job retrieval
      selectInput("prog", "SLiMSuite REST program:", c("zen"),selected=settings$prog),
      textInput("jobid", "REST Server Job ID:", settings$jobid),
      textInput("password", "Job Password [Optional]:", ""),
      textInput("restout", "REST Output to retrieve:", "restkeys"),
      selectInput("restformat", "REST Output format:", c("text","csv","tdt"), "text"),
      #numericInput("pvalue", "P-value cutoff (0-1):", settings$pvalue, min = 0, max = 1, step = 0.01),
      
      actionButton("retrieve", "Retrieve Job"),

      # These parameters alter plot appearance but do not change the data
      checkboxInput("showini", "Show INI file content (Run tab)", value=FALSE),
      checkboxInput("showout", "Show REST output keys (Run tab)", value=TRUE),
      checkboxInput("showdesc", "Show server description", value=FALSE),
      checkboxInput("showinfo", "Show explanation of server output", value=TRUE),
      
      # HTML content to include at end. Contains version number.
      htmlOutput("footer")
      #i# Can include static HTML like this: includeHTML("footer.html")
    ),
    
    # Main output panel
    mainPanel(
      # Overall run info (data$intro)
      verbatimTextOutput("intro"),
      # Tabs of different REST outputs
      tabsetPanel(
        tabPanel("Run", 
                 verbatimTextOutput("status"),
                 conditionalPanel(
                   condition = "input.showini == true",
                   verbatimTextOutput("ini")
                 ),
                 conditionalPanel(
                   condition = "input.showout == true",
                   htmlOutput("restkeys")
                 )
        ),
        tabPanel("Output", 
                 conditionalPanel(
                   condition = "input.restformat == 'text'",
                   #!# Add title then verbatimTextOutput("restout"),
                   htmlOutput("restout")
                 ),
                 conditionalPanel(
                   condition = "input.restformat != 'text'",
                   dataTableOutput("restoutTable")
                 )
        ),
        tabPanel("Wisdoms", 
                 htmlOutput("wisdoms")
        ),
        tabPanel("Log", dataTableOutput("log")),

        tabPanel("Warnings", 
                 verbatimTextOutput("warnings"),
                 textOutput("errortext",container = span),
                 verbatimTextOutput("errors")
        ),

        tabPanel("Help",htmlOutput("help")),
        tabPanel("REST",htmlOutput("outfmt"))
        
        #tabPanel("Plot", plotOutput("simPlot",height=600) ),

        #tabPanel("Summary", verbatimTextOutput("summary")),
        # Main caption based on plotted data
        #tabPanel("Caption", textOutput("caption",container = span)),
        #tabPanel("Table", tableOutput("datatable")),
      ),
      #X#dataTableOutput("sigtable"),
      # Markdown description of plots
      
      #standard = c("status", "version", "ini", "log", "warnings", "error")
      #i# Zen-specific &rest=X outputs
      #noheadtdt = c("wisdoms")
      
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
      
    )
  )
  #<<<# End of ui.R code block #<<<#
))


############### ::: SERVER CODE ::: ##################
#i# server.R defines server logic

