library(shiny)
library(shinybusy)
library(shinythemes)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(shinyalert)
library(shinydashboard)
library(shinyFiles)
library(spatstat.geom)
library(dipsaus)

reslocal <- NULL

pathrdata <- NULL
rdata <- NULL
coords <- NULL
localOpt <- NULL

options(spinner.color="#f06313", spinner.color.background="#ffffff", spinner.size=1)

getIndex <- function(sample){
  
  index <- which(names(reslocal$allTracks.processed)==sample)
  return (index)
  
}

getSamples <- function() {
  
  if( !is.null(reslocal)){
 
    return (names(reslocal$allTracks.processed))}
  
}
#######################################################################

######################### INTERFACE ###################################

#######################################################################

ui <- navbarPage(id="nav_page",
  title="ASCAT.scFit",   theme = shinytheme("united"),
                 tabPanel("Welcome",
                          busy_start_up(
                            loader = tags$img(
                              src = "Loader.gif",
                              width = 400
                            ),
                            text = "Loading ...",
                            color = "#ba4a00",
                            timeout = 3000,
                            background = "white",
                            mode = "auto"
                          ),
                          tags$head(
                            tags$style(HTML(
                              ".hover-effect {color: #FFFFFF ; border-radius: 30px; height:150px; width:300px; background-color: #3e0533; border-color: #6a0144; padding:5px; font-size:4vh; border-width: 5px;  opacity: 0.90;
                              }",
           ".hover-effect:hover {color: #FFFFFF ; border-radius: 30px; background-color: #3e0533; border-color: #6a0144; padding:5px; font-size:5vh; border-width: 5px;
                        opacity: 1;}",
           ".button-container {
             width: 300px;}",
           ".click-effect:active {
      color: #FFFFFF ; border-radius: 30px; background-color: #3e0533; border-color: #6a0144; padding:5px; font-size:5vh; border-width: 5px;
                        opacity: 1;
    
    }",
          "code {
                display:block;
                padding:9.5px;
                margin: auto;
                width: 1100px;
                height: 500px;
                font-size:14px;
                color: #3e0533;
                line-height:4px;
                word-break:break-all;
                word-wrap:break-word;
                white-space:pre-wrap;
                background-color:#FFFFFF;
                border:6px solid #3e0533;
                border-radius:4px; 
            }",
                "pre {
                  display:block;
                  padding:9.5px;
                  margin: auto;
                  width: 350px;
                  height: 100px;
                  font-size:14px;
                  color: #833e03;
                    line-height:5px;
                  word-break:break-all;
                  word-wrap:break-word;
                  white-space:pre-wrap;
                  background-color:#fae6d4;
                    border:4px solid #833e03;
                  border-radius:4px; 
                }",
                "em {
                  display:block;
                  padding:9.5px;
                  margin: auto;
                  width: 800px;
                  height: 95px;
                  font-size:14px;
                  color: #833e03;
                    line-height:5px;
                  word-break:break-all;
                  word-wrap:break-word;
                  white-space:pre-wrap;
                  background-color:#fae6d4;
                    border:4px solid #833e03;
                  border-radius:4px; 
                }"
                                            ))),
                #######################################################################
                
                ######################### WELCOME TAB ################################
                
                #######################################################################
                img(src='WelcomeImage.png', align = "center", style="width:100%; max-width:100%; position: absolute; z-index:-1;"),
                 div(style = "height:70px"),  
                # ba4a00
                h1(strong("ASCAT.sc Ploidy Modifier",style={'color: #5b016a; font-family: arial black ,sans-serif; text-shadow:
  # # 0 0 7px #fff,
  # # 0 0 10px #fff,
  # # 0 0 21px #fff,
  # # 0 0 42px ##5b016a,
  # # 0 0 82px ##5b016a;  font-size: 80px'}), align="center"), br(), br(), br(),
                 fluidRow(column(6, offset = 3, align="center", dropdown(code(h4(strong("1. Choose the data", style={'font-family: Arial'}), align = "left"),
                    p("Please choose the ASCAT.sc rdata object to load. Each sample in the object can be viewed and modified.", align = "left"),  
                    p("Please note that large datasets might take a while to load.", align = "left"),
                   h4(strong("2. Modify the profiles \n ", style={'font-family: Arial; align-text: left'}), align = "left"),
                   p("The purity & ploidy of the whole sample can be changed by clicking on the sunrise plot.", align = "left"), 
                   p("You can additionally shift the sample ploidy, or modify the copy number of 2 different segments.", align = "left"),
                   p("Modifications can be done either by numeric input or directly on the Modifier Station plot.", align = "left"),
                  h4(strong("3. Save the modified profiles ", style={'font-family: Arial'}), align = "left"),
                  p("Once you are happy with the results, you can save all the profiles in text format, or save the entire ASCAT.sc rdata object.", align = "left")),
                                                                         size= "lg", circle= FALSE, status = "info", label = "Help", width = "1200px", inputId="help"))),
                     br(), div(style = "height:50px"),
                fluidRow(column(6, align="center", offset=3, 
                                            shinyFilesButton("get_file", "Load data" ,
                                                             title = "Choose the ASCAT.sc rdata object to load", multiple = FALSE,
                                                             buttonType = "primary", class = NULL, style = "width:60px, height: 50px")
                                            
                     ), br(), br(), br(),div(style = "height:200px"),
                     fluidRow(column(6, align="center", offset = 3,
                                    
                                     actionButtonStyled("start",
                                              label="Start",
                                              class="hover-effect click-effect",
                                              width = "300px"
                                             ))
                                
                                 )
               
                )),
                 
  #######################################################################
  
  ######################### MODIFIER TAB ################################
  
  #######################################################################
  
  tabPanel("Modifier",  shinyWidgets::useShinydashboard(), 
                fluidRow(box(width=12, title="Original", status="warning", solidHeader=TRUE, 
                              column(width=8,withSpinner(plotOutput("profile"),type=3)),
                              column(width=4,plotOutput("sunrise1")))
                  
                 ),
                 useShinyjs(),
                 useShinyalert(),
           fluidRow(column(width=2, offset = 1,
                   dropdownButton(
                     tags$h3("Choose Sample"),
                     br(),
                  selectInput(
                   "samples",
                   label= NULL,
                   choices= getSamples(),
                   selected = NULL,
                   multiple = FALSE,
                   selectize = FALSE
                 ), size= "lg", circle= FALSE, status = "info", icon = icon("list",verify_fa = FALSE), width = "450px",
                label= "Select sample",
                 actionButton("view", label = "View", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 2px")
                 
                   )
                 ),
                
                column(width = 2, offset = 1,
                       dropdownButton(tags$h3("Choose ploidy and purity on the sunrise graph"),
                                      h4("To change the ploidy and purity of the whole sample directly on the Working station profile, please click the point on the sunrise plot corresponding to the desired ploidy/purity values. When 'Automatic optima' is enabled; the closest optimal purity/ploidy pair will be selected."),
                                      checkboxInput("optima", "Automatic optima", TRUE),
                                       size= "lg", circle= FALSE, status = "info", label = "Modify purity & ploidy", width = "450px"
                       )
                       
                ),
                
                column(width = 5,  
                       dropdown(tags$h3("Choose a different way to modify the profile:"),
                                      
                                dropdownButton(tags$h3("Modify the copy number of 2 segments"),br(),
                                               selectInput(
                                                 "Chr1",
                                                 "Choose first chromosome",
                                                 c(1:22, "X", "Y"),
                                                 selected = NULL,
                                                 multiple = FALSE,
                                                 selectize = FALSE,
                                                 width = "75%"
                                                 
                                               ),
                                               textInput("cn1", "Choose first copy number", value = "", placeholder = NULL, width = "75%"),
                                               selectInput(
                                                 "Chr2",
                                                 "Choose second chromosome",
                                                 c(1:22, "X", "Y"),
                                                 selected = NULL,
                                                 multiple = FALSE,
                                                 selectize = FALSE,
                                                 width = "75%"
                                               ),
                                               textInput("cn2", "Choose second copy number", value = "",  placeholder = NULL, width = "75%"),
                                               actionButton("modify", label = "Apply", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                                               size= "lg", circle= FALSE, status = "info",right= TRUE, label = "Refit segments", width = "450px"),
                                br(),
                                dropdownButton(tags$h3("Modify segment on graph"),
                                               h4("To modify the copy number of 2 different segments directly on the Working station profile, please click on the segment you wish to modify, then on its desired position. Repeat for the second segment, then click on 'Apply'"),
                                               actionButton("refit", label = "Apply", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                                               size= "lg", circle= FALSE, status = "info", right= TRUE, label = "Refit segments on graph", width = "450px"
                                ),
                                br(),
                                      dropdownButton(
                                        sliderInput("ploidy", "Shift ploidy by:", -3, 4, 1, step = 1, round = FALSE,
                                                    ticks = TRUE, animate = FALSE,
                                                    width = NULL, sep = ",", pre = NULL, post = NULL),
                                        actionButton("shift", label = "Apply", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                                        size= "lg", circle= FALSE, status = "info", right= TRUE,label = "Shift sample ploidy", width = "450px"
                                      ),
                                br(),
                                      dropdownButton(tags$h3("Shift ploidy on graph"),
                                                     h4("To shift the ploidy of the whole sample directly on the Working station profile, please click on a point with the y axis position corresponding to the desired ploidy value, then click on 'Apply'"),
                                                     actionButton("shift_graph", label = "Apply", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                                                     size= "lg", circle= FALSE, status = "info", right= TRUE,label = "Shift ploidy on graph", width = "450px",up=TRUE
                                      ),
                                      
                                      
                                      size= "lg", circle= FALSE, status = "info", label = "Additional tools", width = "800px", inputId="menu", up=TRUE))
                
                ), br(),
            
    
                fluidRow(box(width=12, title="Modifier Working Station", status="warning", solidHeader=TRUE, column(width=8, withSpinner(plotOutput("profile2", click="profile2_click"),type=3)), column(width = 4, plotOutput("sunrise2", click = "sunrise2_click")))),
                 fluidRow(column(width = 3,  offset=1, actionButton("discard", label = "Reset profile",  icon = icon("arrows-rotate",verify_fa = FALSE), style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px")),
                          column(width = 3, offset=1,  downloadButton("savetxt", label = "Save profiles", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px")),
                                                    column(width = 3, offset=1, downloadButton("save", label = "Save .Rda", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"))), br(), br()))
 
#######################################################################

######################### SERVER ######################################

#######################################################################


server <- function(input, output, session) {
  
  vals <- reactiveVal()
  volumes = getVolumes()
  coords <- reactiveValues(x=NULL,y=NULL)
  output$profile <- renderPlot(NULL)
  output$profile2 <- renderPlot(NULL)
  
  observeEvent(input$profile2_click, {
    coords$x <- c(coords$x,input$profile2_click$x)
    coords$y <- c(coords$y,input$profile2_click$y)
  })
  
  observeEvent(input$sunrise2_click, {
    coords$x <- input$sunrise2_click$x
    coords$y <- input$sunrise2_click$y
  })
  
  sampleName <- reactive({
    input$samples
  })
  optValue <- reactive({input$optima})
  result <- reactive({ reslocal })
  
  shiftv <- reactive({
    input$ploidy
  })
  
  chrs <- reactive({
    list(input$Chr1,input$Chr2, input$cn1, input$cn2)
  })
 
  #########################################################################
  
  ######################## FILE CHOOSER ###################################
  
  #########################################################################
  
  observe({
    shinyFileChoose(input, "get_file", roots = volumes, session = session)

    if(!is.null(input$get_file)){
     
      file_selected<-parseFilePaths(volumes, input$get_file)
      pathrdata <<- file_selected
      
     
    }
  })
  
  observe({ toggle(id="start", condition=(input$get_file>=1))})
  
  #########################################################################
  
  ######################## START ##########################################
  
  #########################################################################
  
  observeEvent(input$start, {
     updateNavbarPage(session=session,
                     inputId="nav_page",
                     selected="Modifier")
    tryCatch({
      
      
      filepath <<- as.character(pathrdata$datapath)
      
      load(filepath)
      reslocal <<- res
     
      updateSelectInput(session, "samples", label=NULL,
                        choices=getSamples())
      if(!"allSolutions.refitted.manual"%in%names(reslocal))
      {
        reslocal$allProfiles.refitted.manual <<- reslocal$allProfiles.refitted.auto
        reslocal$allSolutions.refitted.manual <<- reslocal$allSolutions.refitted.auto
       
      }
      output$profile <- renderImage({
        
        outfile <- tempfile(fileext='.png')
        
        
        png(outfile, width=950, height=400)
        
        
        
        plotSolution(reslocal$allTracks.processed[[1]],
                       purity=reslocal$allSolutions.refitted.auto[[1]]$purity,
                       ploidy=reslocal$allSolutions.refitted.auto[[1]]$ploidy,
                       ismale=if(!is.null(reslocal$sex)) reslocal$sex[[1]]=="male" else "female",
                       gamma=.55,
                       sol=reslocal$allSolutions[[1]])
        
        dev.off()
        
        
        list(src = outfile,
             alt = "Original profile")
      }, deleteFile = TRUE)
     
      output$profile2 <- renderPlot({isolate(plot2 <- plotSolution(reslocal$allTracks.processed[[1]],
                                                                   purity=reslocal$allSolutions.refitted.manual[[1]]$purity,
                                                                   ploidy=reslocal$allSolutions.refitted.manual[[1]]$ploidy,
                                                                   ismale=if(!is.null(reslocal$sex)) reslocal$sex[[1]]=="male" else "female",
                                                                   gamma=.55,
                                                                   sol=reslocal$allSolutions[[1]]))
      })
     
      output$sunrise1 <- renderPlot({isolate(plotSunrise(reslocal$allSolutions.refitted.auto[[1]]))})
      output$sunrise2 <- renderPlot({isolate( localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[1]], localMinima=TRUE))}) 
     
      
      removeModal()
      },
    error=function(e) {
     
    })
   
  })
  
  #########################################################################
  
  ######################## DISCARD/KEEP ###################################
  
  #########################################################################
  
  observeEvent(input$discard, {
    
    vals(input$samples)
    
    if(! is.null(sampleName())){
      index <- getIndex(sampleName())
     
      output$profile2 <- renderPlot({isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                                          purity=reslocal$allSolutions.refitted.auto[[index]]$purity,
                                                          ploidy=reslocal$allSolutions.refitted.auto[[index]]$ploidy,
                                                          ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                                          gamma=.55,
                                                          sol=reslocal$allSolutions[[index]]))})
      
      output$sunrise2 <- renderPlot({isolate(plotSunrise(reslocal$allSolutions.refitted.auto[[index]], localMinima=TRUE))})
      reslocal$allProfiles.refitted.manual[[index]] <<-  reslocal$allProfiles.refitted.auto[[index]] 
      reslocal$allSolutions.refitted.manual[[index]] <<-  reslocal$allSolutions.refitted.auto[[index]] 
      
    }
    else{
      return (NULL)
    }
    
  })
  
  
  #########################################################################
  
  ######################## DOWNLOAD #######################################
  
  #########################################################################
  
  
  output$save <- downloadHandler(
    filename = function() {
      "result_manualfitting.Rda"
    },
    content = function(file) {
     
      showModal(modalDialog(div(tags$b("Loading...", style = "color: steelblue;")), footer=NULL))
      on.exit(removeModal())
      
      res <- reslocal
      save(res, file=file)
      
    }
  )
  output$savetxt <- downloadHandler(
    filename = function(){
      "Profiles_txt.zip"
    },
    content = function(file){
      
      showModal(modalDialog(div(tags$b("Loading...", style = "color: steelblue;")), footer=NULL))
      on.exit(removeModal())
      
      temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_directory)
      
      for ( i in 1:length(reslocal$allProfiles.refitted.manual)){
       
       
        
        write.table(reslocal$allProfiles.refitted.manual[[i]],
                    quote=F,
                    sep="\t",
                    col.names=T,
                    row.names=F,
                    file=paste0(temp_directory,"/",
                                paste0(names(reslocal$allTracks.processed)[i],"-manual_refit"),
                                ".ASCAT.scprofile.txt"))
      }
      
      zip::zip(
        zipfile = file,
        files = dir(temp_directory),
        root = temp_directory
      )
      
    },
    contentType = "application/zip"
    
  )
  
  
  #########################################################################
  
  ######################## SUNRISE ########################################

  #########################################################################
  
 
  observeEvent(input$sunrise2_click,{
    vals <- as.numeric(chrs())
    
    
    if(! is.null(chrs())){
      index <- getIndex(sampleName())
      
      purity <- NULL
      ploidy <- NULL
      
      
      tryCatch(
        {
          
          solution <- reslocal$allSolutions[[1]]
          errs <- solution$errs
          errs <- errs-min(errs)
          errs.max <- max(solution$errs[!is.infinite(solution$errs)])
          errs[is.infinite(errs)] <- errs.max
          errs <- errs/errs.max
          errs <- errs[rev(seq_len(nrow(errs))), ]
          
          purity <- coords$y
          ploidy <- coords$x
         
          if (optValue()){
            best <- 1
            dist <- crossdist(ploidy, purity, localOpt$bao[best, 2]/ncol(errs), 1-localOpt$bao[best, 1]/nrow(errs))
            
            for (i in 1:nrow(localOpt$bao)){
              
              dist2 <- crossdist(ploidy, purity, localOpt$bao[i, 2]/ncol(errs), 1-localOpt$bao[i, 1]/nrow(errs))
              if (dist >= dist2){
                
                dist <- dist2
                best <- i
              }
              
            }
            
            ploidy <- localOpt$bao[best, 2]/ncol(errs)
            purity <- 1-localOpt$bao[best, 1]/nrow(errs)
            
          }
          
          ploidy <- as.numeric(colnames(errs)[as.numeric(ploidy)*ncol(errs)])
          purity <- as.numeric(rownames(errs)[(1-as.numeric(purity))*nrow(errs)])
         
          
          reslocal$allProfiles.refitted.manual[[index]] <<- getProfile(fitProfile(tracksSingle = reslocal$allTracks.processed[[index]], purity, ploidy, ismale=reslocal$sex[index]=="male"))
         
         
          reslocal$allSolutions.refitted.manual[[index]]$ploidy <<- ploidy
         
          reslocal$allSolutions.refitted.manual[[index]]$purity <<- purity
          
          output$profile2 <- renderPlot({isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                                              purity=purity,
                                                              ploidy=ploidy,
                                                              gamma=.55))})
        
          output$sunrise2 <- renderPlot({isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))})
          coords$y <<- NULL
          coords$x <<- NULL
             },
        error=function(e) {
          message('An Error Occurred')
          print(e)
          shinyalert("Error", "Cannot fit profile: ploidy<0 or purity ∉ [0,1]. Please choose different values", type = "error")
        },
        warning=function(w) {
          message('A Warning Occurred')
          print(w)
          shinyalert("Warning", "New solution is ambiguous: reverted to old one", type = "error")
        })
    }
    else{
      return (NULL)
    }
    
  })
  #########################################################################
  
  ######################## MODIFY SEGMENTS ON GRAPH #######################
  
  #########################################################################
  
  observeEvent(input$refit,{
    vals <- as.numeric(chrs())
    
    
    if(! is.null(chrs())){
      index <- getIndex(sampleName())
      input$profile2_click
      chr1 <- NULL
      chr2 <- NULL
      breaks <- c(0, cumsum(sapply(reslocal$allTracks.processed[[1]]$lSegs, function(x) max(x$output$loc.end))/1e+06))
     
      tryCatch(
        {
          for (i in 1:length(breaks)){
            if (coords$x[1] < breaks[i]) {
              chr1 <- i-1
              break} 
           
          }
          for (i in 1:length(breaks)){
            
            if (coords$x[3] < breaks[i]) {
              chr2 <- i-1 
              break} 
          }
         
          y2 <- round(coords$y[2], digits=0)
          y4 <- round(coords$y[4], digits=0)
          
          rescopy <- reslocal
          rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles.refitted.auto <- rescopy$allProfiles.refitted.manual
          
          rescopy <- run_any_refitProfile(rescopy,
                                          sample_indice=index,
                                          chr1=chr1,
                                          ind1=NA,
                                          n1=y2,
                                          chr2=chr2,
                                          ind2=NA,
                                          n2=y4,
                                          CHRS=c(1:22,"X","Y"),
                                          outdir="./www",
                                          gridpur=seq(-.05,.05,.01),
                                          gridpl=seq(-.1,.2,.01))
          rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
          rescopy$allProfiles.refitted.auto <- reslocal$allProfiles.refitted.auto
          
          reslocal <<- rescopy
          output$profile2 <- renderPlot({isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                                              purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                                              ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                                              ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                                              gamma=.55,
                                                              sol=reslocal$allSolutions[[index]]))})
          
          output$sunrise2 <- renderPlot({isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))})
          coords$y <<- NULL
          coords$x <<- NULL
          },
        error=function(e) {
          message('An Error Occurred')
          print(e)
          shinyalert("Error", "Cannot fit profile: ploidy<0 or purity ∉ [0,1]. Please choose different values", type = "error")
        },
        warning=function(w) {
          message('A Warning Occurred')
          print(w)
          shinyalert("Warning", "New solution is ambiguous: reverted to old one", type = "error")
        })
    }
    else{
      return (NULL)
    }
    
  })
  #########################################################################
  
  ######################## MODIFY SEGMENTS ################################
  
  #########################################################################
  
  observeEvent(input$modify,{
    vals <- as.numeric(chrs())
    
    
    if(! is.null(chrs())){
      index <- getIndex(sampleName())
      tryCatch(
      {
        rescopy <- reslocal
        rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
        rescopy$allProfiles.refitted.auto <- rescopy$allProfiles.refitted.manual
        
        rescopy <- run_any_refitProfile(rescopy,
                                  sample_indice=index,
                                  chr1=vals[[1]],
                                  ind1=NA,
                                  n1=vals[[3]],
                                  chr2=vals[[2]],
                                  ind2=NA,
                                  n2=vals[[4]],
                                  CHRS=c(1:22,"X","Y"),
                                  outdir="./www",
                                  gridpur=seq(-.05,.05,.01),
                                  gridpl=seq(-.1,.2,.01))
        rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
        rescopy$allProfiles.refitted.auto <- reslocal$allProfiles.refitted.auto
        
        reslocal <<- rescopy
      output$profile2 <- renderPlot({isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                                          purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                                          ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                                          ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                                          gamma=.55,
                                                          sol=reslocal$allSolutions[[index]]))})
      
      output$sunrise2 <- renderPlot({isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))})
     
        },
      error=function(e) {
        message('An Error Occurred')
        print(e)
        shinyalert("Error", "Cannot fit profile: ploidy<0 or purity ∉ [0,1]. Please choose different values", type = "error")
      },
      warning=function(w) {
        message('A Warning Occurred')
        print(w)
        shinyalert("Warning", "New solution is ambiguous: reverted to old one", type = "error")
      })
    }
    else{
      return (NULL)
    }
    
  })
  #########################################################################
  
  ######################## SHIFT ON GRAPH #################################
  
  #########################################################################
  
  observeEvent(input$shift_graph,{
    vals <- as.numeric(chrs())
    input$profile2_click
    
    if(! is.null(chrs())){
      
      
      tryCatch(
        {
          index <- getIndex(sampleName())
          
          ypos <- as.numeric(round(coords$y[1], digits=0))
          ploidy <- round(as.numeric(reslocal$allSolutions.refitted.manual[[index]]$ploidy), digits=0)
          
          shiftp <- ypos - ploidy
          rescopy <- reslocal
          rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles.refitted.auto <- rescopy$allProfiles.refitted.manual
          
          rescopy <- run_any_refitProfile_shift(rescopy,
                                                sample_indice=index,
                                                shift=shiftp,
                                                CHRS=c(1:22,"X","Y"),
                                                outdir="./www",
                                                gridpur=seq(-.05,.05,.01),
                                                gridpl=seq(-.1,.2,.01))
          
          rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
          rescopy$allProfiles.refitted.auto <- reslocal$allProfiles.refitted.auto
          
          reslocal <<- rescopy
          
          output$profile2 <- renderPlot({isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                                              purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                                              ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                                              ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                                              gamma=.55,
                                                              sol=reslocal$allSolutions[[index]]))})
          
          output$sunrise2 <- renderPlot({isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))})
          coords$y <<- NULL
          coords$x <<- NULL
          
        },
        error=function(e) {
          message('An Error Occurred')
          print(e)
          shinyalert("Error", "Cannot fit profile. Please choose different values", type = "error")
        },
        warning=function(w) {
          message('A Warning Occurred')
          print(w)
          shinyalert("Warning", "New solution is ambiguous: reverted to old one", type = "error")
        }
      )
      
    }
    else{
      return (NULL)
    }
    
  })
  #########################################################################
  
  ######################## SHIFT ##########################################
  
  #########################################################################
  
  observeEvent(input$shift,{
    vals <- as.numeric(chrs())
    
   
    if(! is.null(chrs())){
      index <- getIndex(sampleName())
      
      
      withCallingHandlers(
        {
          rescopy <- reslocal
          rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles.refitted.auto <- rescopy$allProfiles.refitted.manual
          
          shiftp <- as.numeric(shiftv())
          
        
            rescopy <- run_any_refitProfile_shift(rescopy,
                                                  sample_indice=index,
                                                  shift=shiftp,
                                                  CHRS=c(1:22,"X","Y"),
                                                  outdir="./www",
                                                  gridpur=seq(-.05,.05,.01),
                                                  gridpl=seq(-.1,.2,.01))
            
            rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
            rescopy$allProfiles.refitted.auto <- reslocal$allProfiles.refitted.auto
            
            reslocal <<- rescopy
            
            output$profile2 <- renderPlot({plotSolution(reslocal$allTracks.processed[[index]],
                                                                purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                                                ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                                                ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                                                gamma=.55,
                                                                sol=reslocal$allSolutions[[index]])})
            
            output$sunrise2 <- renderPlot({plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE)})
          
        },
        error=function(e) {
           shinyalert("Error", "Cannot fit profile. Please choose different values", type = "error")
        },
        warning=function(w) {
         
          shinyalert("Warning", "Cannot fit profile. Please choose different values", type = "error")
        }
      )
     
   
    }
    else{
      return (NULL)
    }
    
  })
  #########################################################################
  
  ######################## VIEW ###########################################
  
  #########################################################################
  
  observeEvent(input$view,{
    vals(input$samples)
    
    if(! is.null(sampleName())){
      index <- getIndex(sampleName())
    
              
      output$profile <- renderImage({
        
        outfile <- tempfile(fileext='.png')
        
        
        png(outfile, width=950, height=400)
        
        
        
        plotSolution(reslocal$allTracks.processed[[1]],
                       purity=reslocal$allSolutions.refitted.auto[[1]]$purity,
                       ploidy=reslocal$allSolutions.refitted.auto[[1]]$ploidy,
                       ismale=if(!is.null(reslocal$sex)) reslocal$sex[[1]]=="male" else "female",
                       gamma=.55,
                       sol=reslocal$allSolutions[[1]])
        
        dev.off()
        
        
        list(src = outfile,
             alt = "Original profile")
      }, deleteFile = TRUE)
       output$profile2 <- renderPlot({isolate(plot2 <- plotSolution(tracksSingle=reslocal$allTracks.processed[[index]],
                                                                    purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                                                    ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                                                    ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                                                    gamma=.55,
                                                                    sol=reslocal$allSolutions[[index]]))
       })
       output$sunrise1 <- renderPlot({isolate(plotSunrise(reslocal$allSolutions.refitted.auto[[index]]))})
       output$sunrise2 <- renderPlot({isolate( localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))})
       
       localOpt <<- localOpt$bao
     
  }
    else{
      return (NULL)
    }
    })
  
} 

shinyApp(ui = ui, server = server)
