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
      #i# This program selection is used experimentally to control output. It might be possible to remove.
      #i# Would be better to change to select an output style
      selectInput("prog", "SLiMSuite REST program:", c("None"),selected="None"),
      # Set parameters for REST job retrieval
      textInput("jobid", "REST Server Job ID:", settings$jobid),
      textInput("password", "Job Password [Optional]:", ""),
      selectInput("restout", "REST Output to retrieve:", c("status"), "status"),
      #X#textInput("restout", "REST Output to retrieve:", "status"),
      #!# Change to plain, text and table
      #selectInput("restformat", "REST Output format:", c("none","text","table","plot","class"), "none"),
      selectInput("restformat", "REST Output format:", c("none","text","table","plot"), "none"),
      
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
      conditionalPanel(
        condition = "input.prog != 'None'",
        htmlOutput("resultsChoice")
      ),
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
        #i# Output panel for formatted REST results output
        tabPanel("Results", 
                 #i# none = Plain text format
                 conditionalPanel(
                   condition = "input.restformat == 'none'",
                   htmlOutput("restoutPlain")
                 ),
                 #i# text = Regular text format
                 conditionalPanel(
                   condition = "input.restformat == 'text'",
                   #!# Add title then verbatimTextOutput("restout"),
                   htmlOutput("restoutText")
                 ),
                 #i# table = Data Table format
                 conditionalPanel(
                   condition = "input.restformat == 'table'",
                   #!# Add title then verbatimTextOutput("restout"),
                   dataTableOutput("restoutTable")
                 ),
                 #i# plot = Server-specific plot
                 conditionalPanel(
                   condition = "input.restformat == 'plot'",
                   #!# Add title then verbatimTextOutput("restout"),
                   htmlOutput("restoutTitle")
                   #!# htmlOutput("restoutPlot")
                 ),
                 #i# class = Server output classes
                 conditionalPanel(
                   condition = "input.restformat == 'class'",
                   #!# Add title then verbatimTextOutput("restout"),
                   htmlOutput("restoutClass")
                   #!# htmlOutput("restoutPlot")
                 )
        ),
        #X# Removing Wisdoms in favour of dynamic Output
        #tabPanel("Wisdoms", 
        #       htmlOutput("wisdoms")
        #),
        tabPanel("Log", dataTableOutput("log")),

        tabPanel("Warnings", 
                 verbatimTextOutput("warnings"),
                 textOutput("errortext",container = span),
                 verbatimTextOutput("errors")
        ),

        tabPanel("Help",htmlOutput("help")),
        tabPanel("REST",htmlOutput("outfmt")),
        tabPanel("Retrieve",htmlOutput("retrieve"))
        
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

