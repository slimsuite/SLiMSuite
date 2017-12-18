#i# See the `main.R` file for version, contact and license information.
#i# The `main.R` script loads libraries and contains the initial parameter settings and functions.
source("main.R")

# This is the code that goes inside the server object
shinyServer(function(input, output, session) {

  ### SECTION 1 - Setup a reactive data list for the server data
  adata <- reactiveValues(
      data = setupData()
  )
  #i# Can update in a render function using adata$data = ...

  
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
          updateSelectInput(session, "prog",
                            label = "SLiMSuite REST program::",
                            choices = c(adata$data$prog,"None"),
                            selected = adata$data$prog
          )
          rkeys = c()
          for(rkey in names(adata$data)){
            if(!rkey %in% settings$stdkeys){
              rkeys = c(rkeys,rkey)
            }
          }
          updateSelectInput(session, "restout",
              label = "REST Output to retrieve:",
              #choices = adata$data$restkeys,
              #choices = names(adata$data),
              choices = rkeys,
              selected = "restkeys"
          )
        })
      })
    }
    return(paste(as.character(adata$data$status),sep="\n",collapse="\n"))
  })
  #i# Additional text output reporting what is in the Results tab
  output$resultsChoice <- renderUI({
    if(input$retrieve > 0){
      updateSelectInput(session, "restformat",
                        selected = getRestFormat(input$restout,pure=FALSE)
      )
      myhtml = paste0("<p>Selected output in <b>Results</b> tab: <code>",input$restout,"</code></p>")
      return(HTML(myhtml))
    }
  })
  
  ### SECTION 3 - Output tabs: data rendering
  ### Standard Server Outputs
  # Verbatim text outputs
  output$intro <- renderText({
    #return(paste(as.character(adata$data$status[,1]),sep="\n")) 
    return(paste(as.character(adata$data$intro),sep="\n",collapse="\n"))
  })
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
    updateSelectInput(session, "restformat",
                      selected = getRestFormat(input$restout,pure=FALSE)
    )
    pretext = paste(as.character(adata$data$restkeys),sep="\n",collapse="\n")
    myhtml = c(paste0("<h2>",adata$data$prog," Outputs</h2>"),
               paste0("<p>Selected output in <b>Results</b> tab: <code>",input$restout,"</code></p>"),
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
  
  # HMTL output or URL linking to Job on REST server
  output$retrieve <- renderUI({
    if(isJobID(input$jobid)){
      joburl = paste0(settings$resturl,"retrieve&jobid=",input$jobid,"&password=",input$password,"&rest=parse")
      joblink = paste0("<p>Click here to open the job on the SLiMSuite REST server:</p>\n<p><a href=\"",joburl,"\">",joburl,"</a></p>\n")
      return(HTML(joblink))
    }else{
      return(HTML(paste("<b>ERROR:",input$jobid,"is an invalid JobID.</b>\n"),sep="\n",collapse="\n"))
    }
  })
  
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
  #i# tabPanel("Log", dataTableOutput("log")),
  output$log = renderDataTable({
      adata$data[["log"]]
    },
    rownames=FALSE,
    options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 25)
  )
  
  ### Generic REST server output rendering
  #!# NB. This does not work!
  reactiveOutFmt <- reactive({
    updateSelectInput(session, "restformat",
                      selected = getRestFormat(input$restout,pure=FALSE)
    )
  })
  #i# Plain text output
  output$restoutPlain <- renderUI({
    myhtml = c(paste0("<h3>",input$restout," (",input$restformat,"):</h3>\n"))
    #pretxt = paste(as.character(adata$data[[input$restout]]),sep="\n",collapse="\n")
    pretxt = paste0('<pre>',paste(as.character(adata$data[[input$restout]]),sep="\n",collapse="\n"),'</pre>\n')
    myhtml = c(myhtml,pretxt) 
    HTML(paste0(myhtml),collapse='')
  })
  #i# Regular text output
  output$restoutText <- renderUI({
    #!# Add reading of type from function and dynamically updating the type selection and Output tab label
    myhtml = c(paste0("<h3>",input$restout," (",input$restformat,"):</h3>\n<p>"))
    myhtml = c(myhtml,paste(as.character(adata$data[[input$restout]]),sep="</p>\n<p>",collapse="</p>\n<p>","</p>\n"))
    HTML(paste0(myhtml))
    #return(HTML(paste(as.character(adata$data$wisdoms[,1]),sep="<br/>\n",collapse="<br/>\n")))
  })
  #i# Class type output
  output$restoutClass <- renderUI({
    #!# Add reading of type from function and dynamically updating the type selection and Output tab label
    myhtml = c(paste0("<h3>",input$restout," (",input$restformat,"):</h3>\n<p>"))
    myhtml = c(myhtml,class(adata$data[[input$restout]]),sep="</p>\n<p>",collapse="</p>\n<p>","</p>\n")
    HTML(paste0(myhtml))
    #return(HTML(paste(as.character(adata$data$wisdoms[,1]),sep="<br/>\n",collapse="<br/>\n")))
  })
  #i# Data table output
  output$restoutTable = renderDataTable({
      adata$data[[input$restout]]
    },
    rownames=FALSE,
    options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 25)
  )
  
  
  ### Specifc server output rendering
  #i# Return multi-line HTML coded text
  output$wisdoms <- renderUI({
    myhtml = c("<h2>Random wisdom generator:</h2>\n<p>")
    myhtml = c(myhtml,paste(as.character(adata$data$wisdoms),sep="</p>\n<p>",collapse="</p>\n<p>","</p>\n"))
    HTML(paste0(myhtml))
    #return(HTML(paste(as.character(adata$data$wisdoms[,1]),sep="<br/>\n",collapse="<br/>\n")))
  })
  
  
  
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
  

  #i# Debugging output
  output$debug <- renderUI({
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
