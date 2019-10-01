library(shiny)
library(shinyjs)
library(shinyBS)
library(rintrojs)
library(shinythemes)
library(DT)
library(markdown)
library(SigHotSpotter)
library(visNetwork)

#################################
## UI definition
#################################
analysis_page <- fluidPage(
  theme = shinytheme("spacelab"),
  useShinyjs(),
  # App title
  title='analysis_page',
  h1("Analysis", inline=TRUE),
  shiny::tags$hr(),
  introjsUI(), #Needed for take a tour
  h3(textOutput("header", inline=TRUE),
     div(

       actionButton("tour","Take a tour"),
       actionButton("restart", "Restart"),
       br(),
       br(),
              introBox(
                downloadButton("downloadResults", "Download"),
                data.step = 8,
                data.intro = "Click to download all the generated results"
              ),
       class='rightAlign')
     ),

  div( id ="Sidebar",
       sidebarPanel( #width = 3,
         textOutput("description1", inline = TRUE),
         #TODO: help button
         #bsButton("q1", label = "", icon = icon("question-circle-o", "fa-2x"), style = "default", size = "extra-small"),
         #bsTooltip("q1", 'Naming convention for columns is "Parental...", "Daugher1...", and "Daughter2..." for progenitor and the two daugther cell-types.',
                   #"right", options = list(container = "body")),
         introBox(
         fileInput("cond1_file", "Choose file for condition 1",
                   multiple = FALSE,
                   accept = c("text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv", ".tsv" )),
         data.step = 1,
         data.intro = "Please upload data file corresponding to condition 1"
         ),

         introBox(
         fileInput("cond2_file", "Choose file for condition 2",
                   multiple = FALSE,
                   accept = c("text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv", ".tsv" )),
         data.step = 2,
         data.intro = "Please upload data file corresponding to condition 2"
         ),

         introBox(
         fileInput("de_file", "Choose file for differential expression",
                   multiple = FALSE,
                   accept = c("text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv", ".tsv" )),
         data.step = 3,
         data.intro = "Please upload data file for differential expression of condition 1 compared to condition 2"
         ),
             shiny::tags$hr(),

         introBox(
           selectInput("species", "Species",
                       choices = c("Please select" = "", "Human" = "HUMAN", "Mouse" = "MOUSE") ),
         data.step = 4,
         data.intro = "Select the species. Currently, the supported species are Human and Mouse"
         ),

        br(),
        fluidRow(
          column(6,
         introBox(
                 numericInput("cutoff", "Cutoff:", 50, min = 0, max = 100),
         data.step = 5,
         data.intro = "Select cutoff for percentage of samples in which the gene has to have non-zero expression"
         )
          ),
          column(6,
         introBox(
                 numericInput("pctile", "Percentile:", 90, min = 1, max = 100),
         data.step = 6,
         data.intro = "Select cutoff for top percentile for predicting intermediates"
         )
          )
        ),
        br(),

        fluidRow(
          column(6,
                 introBox(
                   actionButton("button", "Run"),
                   data.step = 7,
                   data.intro = "Click on Run to start the analysis"
                 )
          )
        )
       )
  ),
  # Main panel for displaying data
  mainPanel(

    fluidRow(
      column(5, h3(textOutput('titleText1', inline=TRUE)))
    ),
    fluidRow(
      column(7, class="centerAlign", DT::dataTableOutput("results1_active")),
        #Thisis 1
      column(5, div( id ="div_A1",
        fluidRow(
          visNetworkOutput("network_A1")
        ),
        fluidRow(
          div( id ="divButton_A1", downloadButton("download_A1", "Download network") )
        )
      ), offset = 0)
    ),
    fluidRow(
      column(7, class="centerAlign", DT::dataTableOutput("results1_inactive")),
      column(5, div( id ="div_I1",
        fluidRow(
          visNetworkOutput("network_I1")
        ),
        fluidRow(
          div( id ="divButton_I1", downloadButton("download_I1", "Download network") )
        )
      ), offset = 0)
    ),
    fluidRow(
      column(5, h3(textOutput('titleText2', inline=TRUE)))
    ),
    fluidRow(
      column(7, class="centerAlign", DT::dataTableOutput("results2_active")),
      column(5, div( id ="div_A2",
        fluidRow(
          visNetworkOutput("network_A2")
        ),
        fluidRow(
          div( id ="divButton_A2", downloadButton("download_A2", "Download network") )
        )
      ), offset = 0)
    ),
    fluidRow(
      column(7, class="centerAlign", DT::dataTableOutput("results2_inactive")),
      column(5, div( id ="div_I2",
        fluidRow(
          visNetworkOutput("network_I2")
        ),
        fluidRow(
          div( id ="divButton_I2", downloadButton("download_I2", "Download network") )
        )
      ), offset = 0)
    ),

    width=10
  )

)

help_page <- fluidPage(h1("Help"),title='help_page',navlistPanel(
    tabPanel("Documentation",
             {includeMarkdown('www/Documentation.md')}
    ),
    tabPanel("About",
             {includeMarkdown('www/About.md')}
    ),
    tabPanel("Terms",
             {includeMarkdown('www/Terms.md')}
    ),
    widths = c(2,10))
)
source_page <- fluidPage(h1("Source"),title='source_page',  {includeMarkdown('www/Source.md')} )
documentation_page <- fluidPage(h1("Documentation"),title='documentation_page')
contact_page <- fluidPage(h1("Contact"),title='contact_page')
corner_element <- HTML("<a style=\"color:currentColor;text-decoration:none;\" rel=\"Home\" href=\"http://sighotspotter.lcsb.uni.lu/\">SigHotSpotter</a>")
ui <- shinyUI(
  fluidPage(
    shiny::tags$head( includeCSS('www/default.css'), includeHTML('www/favicon.html')),
    includeHTML('www/header.html'),
    navbarPage(
               corner_element,
               tabPanel("Analysis", analysis_page),

               #TODO: Help and source pages
               tabPanel("Help", help_page),
               tabPanel("Source", source_page),
               #tabPanel("Documentation", documentation_page),
               #tabPanel("Contact", contact_page),
               collapsible = TRUE,
               windowTitle = "SigHotSpotter"
               #TODO: include icon
              # icon = {icon("calendar")}, #{'https://d30y9cdsu7xlg0.cloudfront.net/png/230171-200.png'},
    ),
    shiny::tags$hr(),
    div(class='footer', textOutput("footertext", inline = FALSE))
    ,{includeHTML('www/footer.html')}

  )
)

# Function to plot the Visnetowkr object
# @param visg the visnetwork object
vis.net.plot <- function(visg){
  #hierarchy
  visNetwork(visg$nodes,visg$edges) %>% visNodes(visg, shape="box") %>%
    visIgraphLayout(layout = "layout_as_tree",root="NICHE",flip.y = F) %>%
    visEdges(arrows = "to") %>%  visOptions(highlightNearest = list(enabled =TRUE, degree = 1, hover = T), nodesIdSelection = TRUE)  %>%
    visEdges(smooth = T) %>% visGroups(visg, groupname="int", shape="circle",color="blue") %>%
	visNodes(visg,color="grey") %>%
    visGroups(visg, groupname="upregulated", color = "red",shape="triangle") %>%
    visGroups(visg, groupname="downregulated", color = "green", shape="triangle") %>%
    visPhysics(stabilization = FALSE) %>% visEdges(smooth = FALSE) %>%
    visExport(float = "left")
}

#################################
## Function to update button state based on the status of uploaded files
#################################
.updateButtons <- function(rv) {
  if(rv$species == "" || is.null(rv$cond1_file) || is.null(rv$cond2_file) || is.null(rv$de_file)){
    shinyjs::disable("button")
  } else {
    shinyjs::enable("button")
  }

  if(length(g_results) == 0) {
    shinyjs::disable("downloadResults")
  } else {
    shinyjs::enable("downloadResults")
  }
}

#################################
## Server logic
#################################
server <- function(input, output, session) {

  # Option to increase accepted file size (500Mb)
  options(shiny.maxRequestSize=500*1024^2)

  # Global variables
  g_results <<- list()

  # Create reactive values
  rv <- reactiveValues( species = "", cond1_file = NULL,  cond2_file = NULL, de_file = NULL )

  shinyjs::disable("button")

  output$header <- renderText("Input data")
  #output$description1 <- renderText('Please upload data files containing for the two conditions and the differential analysis')

  output$footertext<- renderText(paste0('SigHotSpotter v',  packageVersion("SigHotSpotter") ) )

  shinyjs::hide(id = "div_A1")

  #################################
  # Observer definitions
  #################################

  #TODO: do we need examples?
  # Example chooser
 if (FALSE) {
  observeEvent(input$example, {
    if(!is.null(input$example))
    {
      shinyjs::reset("dataFile")
      myfiles = strsplit(input$example, ";")[[1]]
      rv$dataFile = system.file("extdata/examples", myfiles[1], package="seesawData")
      updateSelectInput(session, "hairball", selected = myfiles[2])
      .updateButtons(rv)
    }
  } )

 }

  #Species select input
  observeEvent(input$species, {
    rv$species = input$species
    .updateButtons(rv)
  })

  # import data
  observeEvent(input$cond1_file, {
    rv$cond1_file=input$cond1_file
    #updateSelectizeInput(session, "example", selected = "")
    .updateButtons(rv)
  })

  observeEvent(input$cond2_file, {
    rv$cond2_file=input$cond2_file
    #updateSelectizeInput(session, "example", selected = "")
    .updateButtons(rv)
  })

  observeEvent(input$de_file, {
    rv$de_file=input$de_file
    #updateSelectizeInput(session, "example", selected = "")
    .updateButtons(rv)
  })



  # restart button
  observeEvent(input$restart, {
    shinyjs::reset("cond1_file")
    shinyjs::reset("cond2_file")
    shinyjs::reset("de_file")

    #NOTE: reset select input: this does not seem to work locally, but it works on the webserver
    updateSelectInput(session, "species", choices = c("Please select" = "", "Human" = "HUMAN", "Mouse" = "MOUSE") ,
                      selected = character(0))
    #Nor does this
    #shinyjs::reset("species")

    session$reload()
  })

  # Take a tour
  observeEvent(input$tour,{
    introjs(session )
  })

  # Download results

  output$downloadResults <- downloadHandler(
    filename = function() {
      paste("Results-", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file){
      library('openxlsx')

      ## Create a blank workbook
      wb <- createWorkbook()

      # General information
      l_info = list()
      l_info$Software <-  paste0('SigHotSpotter v',  packageVersion("SigHotSpotter") )
      l_info$Condition1 <- rv$cond1_file$name
      l_info$Condition2 <- rv$cond2_file$name
      l_info$Cutoff <- input$cutoff
      l_info$Percentile <- input$pctile

      # Add infosheet to wb
      addWorksheet(wb, 'General')
      writeData(wb, 'General', as.data.frame(unlist(l_info)), rowNames = TRUE, colNames = FALSE)

     # Condition 1
      addWorksheet(wb, 'Condition1')
      condition <- g_results[[1]]$final_score
      condition$Steady_state <- NULL
      colnames(condition) <- c('Signaling hotspots', 'Compatibility score*')
      condition = rbind(condition, c('', ''))
      condition = rbind(condition, c('', '*: Values larger than 0.5 denote active, smaller than 0.5 denote inactive signaling hotspots'))
      writeData(wb, 'Condition1', condition, rowNames = FALSE, colNames = TRUE)

     # Condition 2
      addWorksheet(wb, 'Condition2')
      condition <- g_results[[2]]$final_score
      condition$Steady_state <- NULL
      colnames(condition) <- c('Signaling hotspots', 'Compatibility score*')
      condition = rbind(condition, c('', ''))
      condition = rbind(condition, c('', '*: Values larger than 0.5 denote active, smaller than 0.5 denote inactive signaling hotspots'))
      writeData(wb, 'Condition2', condition, rowNames = FALSE, colNames = TRUE)

      ## Save workbook to working directory
      saveWorkbook(wb, file = file, overwrite = TRUE)
    }
  )

   #Thisis 3
  # Download networks
  # active 1
  output$download_A1<- downloadHandler(
    filename = function() {
      fname = 'network.sif'
      return (fname)
    },
    content = function(file){
      s = input$results1_active_rows_selected
      #to save network in Sif file for cytoscape
      write.table(g_results[[1]]$vis_net_A[[s]]$edges,file,append = F, row.names = F,quote=F)
    }
  )

  # inactive 1
  output$download_I1<- downloadHandler(
    filename = function() {
      fname = 'network.sif'
      return (fname)
    },
    content = function(file){
      s = input$results1_inactive_rows_selected
      #to save network in Sif file for cytoscape
      write.table(g_results[[1]]$vis_net_I[[s]]$edges,file,append = F, row.names = F,quote=F)
    }
  )

  # active 2
  output$download_A2<- downloadHandler(
    filename = function() {
      fname = 'network.sif'
      return (fname)
    },
    content = function(file){
      s = input$results2_active_rows_selected
      #to save network in Sif file for cytoscape
      write.table(g_results[[2]]$vis_net_A[[s]]$edges,file,append = F, row.names = F,quote=F)
    }
  )

  # inactive 2
  output$download_I2<- downloadHandler(
    filename = function() {
      fname = 'network.sif'
      return (fname)
    },
    content = function(file){
      s = input$results2_inactive_rows_selected
      #to save network in Sif file for cytoscape
      write.table(g_results[[2]]$vis_net_I[[s]]$edges,file,append = F, row.names = F,quote=F)
    }
  )


  #################################
  # Output generation
  #################################

  # Reactive dependency on input$button
  observeEvent(input$button, {

    tryCatch(
      {

        # generating the results
        withProgress(message = 'Please wait.',  {
          withProgress(message = 'Condition 1',  {

            g_results[[1]] <<- SigHotSpotter_pipeline (input$species, rv$cond1_file$datapath, input$cutoff, rv$de_file$datapath, input$pctile, invert_DE = FALSE)
            output$results1_active <- DT::renderDataTable({
              g_results[[1]]$trimmed_score_A
            },server = TRUE, selection = "single")

            output$results1_inactive <- DT::renderDataTable({
              g_results[[1]]$trimmed_score_I
            },server = TRUE, selection = "single")
          })

          output$titleText1 <- renderText(paste("Condition1:",rv$cond1_file$name))

          incProgress(0.4)

          withProgress(message = 'Condition 2',  {
            g_results[[2]] <<- SigHotSpotter_pipeline (input$species, rv$cond2_file$datapath, input$cutoff, rv$de_file$datapath, input$pctile, invert_DE = TRUE)
            output$results2_active <- DT::renderDataTable({
              g_results[[2]]$trimmed_score_A
            },server = TRUE, selection = "single")

            output$results2_inactive <- DT::renderDataTable({
              g_results[[2]]$trimmed_score_I
            },server = TRUE, selection = "single")
          })

          output$titleText2 <- renderText(paste("Condition2:",rv$cond2_file$name))
        })

        .updateButtons(rv)
        shinyjs::hide(id = "Sidebar")
        shinyjs::hide(id = "Input")


        #Thisis 2
        shinyjs::show(id = "div_A1")
        shinyjs::hide(id = "divButton_A1")

        shinyjs::show(id = "div_I1")
        shinyjs::hide(id = "divButton_I1")

        shinyjs::show(id = "div_A2")
        shinyjs::hide(id = "divButton_A2")

        shinyjs::show(id = "div_I2")
        shinyjs::hide(id = "divButton_I2")
      }
      , error = function(e){
        showModal(modalDialog(
          title = "Error", e$message,
          footer = modalButton("Dismiss")
        ))
      }
    )

  })

#Thisis 4
  # Plot network when clicking on a result row
  # active 1
  output$network_A1 <- renderVisNetwork({
    s = input$results1_active_rows_selected
    if (length(s)) {
        shinyjs::show(id = "divButton_A1")
      vis.net.plot(g_results[[1]]$vis_net_A[[s]])
    } else {
        shinyjs::hide(id = "divButton_A1")
    }
  })

  # inactive 1
  output$network_I1 <- renderVisNetwork({
    s = input$results1_inactive_rows_selected
    if (length(s)) {
        shinyjs::show(id = "divButton_I1")
      vis.net.plot(g_results[[1]]$vis_net_I[[s]])
    } else {
        shinyjs::hide(id = "divButton_I1")
    }
  })

  # active 2
  output$network_A2 <- renderVisNetwork({
    s = input$results2_active_rows_selected
    if (length(s)) {
        shinyjs::show(id = "divButton_A2")
      vis.net.plot(g_results[[2]]$vis_net_A[[s]])
    } else {
        shinyjs::hide(id = "divButton_A2")
    }
  })

  # inactive 2
  output$network_I2 <- renderVisNetwork({
    s = input$results2_inactive_rows_selected
    if (length(s)) {
        shinyjs::show(id = "divButton_I2")
      vis.net.plot(g_results[[2]]$vis_net_I[[s]])
    } else {
        shinyjs::hide(id = "divButton_I2")
    }
  })

}

#################################
## Running the app (only in interactive environment)
#################################
#if (interactive()) {
  shinyApp(ui, server)
#}
