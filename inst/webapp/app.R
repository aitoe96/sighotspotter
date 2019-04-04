library(shiny)
library(shinyjs)
library(shinyBS)
library(rintrojs)
library(shinythemes)
library(DT)
library(markdown)
library(NicheSIG)

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
         introBox(
        actionButton("button", "Run"),
         data.step = 7,
         data.intro = "Click on Run to start the analysis"
         )
       )
  ),
  # Main panel for displaying data
  mainPanel(

 #TODO: obsoleted, delete if not needed
 #   div( id="Input",
 #        DT::dataTableOutput("contents", width = 900)
 #   ),

    # Panel for displaying results
    # verbatimTextOutput('result_summary'),
    #fluidRow(
      column(5, h3(textOutput('titleText1', inline=TRUE))),
      column(9, class="centerAlign", DT::dataTableOutput("results1_active")),
      column(9, class="centerAlign", DT::dataTableOutput("results1_inactive")),
      column(5, h3(textOutput('titleText2', inline=TRUE))),
      column(9, class="centerAlign", DT::dataTableOutput("results2_active")),
      column(9, class="centerAlign", DT::dataTableOutput("results2_inactive"))
    #)

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
corner_element <- HTML("<a style=\"color:currentColor;text-decoration:none;\" rel=\"Home\" href=\"http://nichesig.lcsb.uni.lu/\">NicheSIG</a>")
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
               windowTitle = "NicheSIG"
               #TODO: include icon
              # icon = {icon("calendar")}, #{'https://d30y9cdsu7xlg0.cloudfront.net/png/230171-200.png'},
    ),
    shiny::tags$hr(),
    div(class='footer', textOutput("footertext", inline = FALSE))
    ,{includeHTML('www/footer.html')}

  )
)

#################################
## Function to update button state based on the status of uploaded files
#################################
.updateButton <- function(rv) {
  if(rv$species == "" || is.null(rv$cond1_file) || is.null(rv$cond2_file) || is.null(rv$de_file)){
    shinyjs::disable("button")
  } else {
    shinyjs::enable("button")
  }
}

#################################
## Function to trim results for display
#################################
.trimResults <- function(results, active = TRUE) {

  res_trimmed = results[,1:2]
  if (active){
    res_trimmed <- head(res_trimmed, 10L)
    colnames(res_trimmed) <- c('Active signaling hotspots', 'Compatibility score')
  } else
  {
    res_trimmed <- tail(res_trimmed, 10L)
    colnames(res_trimmed) <- c('Inactive signaling hotspots', 'Compatibility score')
    res_trimmed <- res_trimmed[order(res_trimmed$'Compatibility score'),]
  }
  res_trimmed[,2] = round(res_trimmed[,2],4)
  rownames(res_trimmed) <- NULL
  return(res_trimmed)
}

#################################
## Server logic
#################################
server <- function(input, output, session) {

  # Option to increase accepted file size (500Mb)
  options(shiny.maxRequestSize=500*1024^2)

  # Global variables
  # seesaw_res <<- NULL
  # TODO: do we need any?


  # Create reactive values
  rv <- reactiveValues( species = "", cond1_file = NULL,  cond2_file = NULL, de_file = NULL )

  shinyjs::disable("button")

  output$header <- renderText("Input data")
  #output$description1 <- renderText('Please upload data files containing for the two conditions and the differential analysis')

  output$footertext<- renderText(paste0('NicheSIG v',  packageVersion("NicheSIG") ) )

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
      .updateButton(rv)
    }
  } )

 }

  #Species select input
  observeEvent(input$species, {
    rv$species = input$species
    .updateButton(rv)
  })

  # import data
  observeEvent(input$cond1_file, {
    rv$cond1_file=input$cond1_file
    #updateSelectizeInput(session, "example", selected = "")
    .updateButton(rv)
  })

  observeEvent(input$cond2_file, {
    rv$cond2_file=input$cond2_file
    #updateSelectizeInput(session, "example", selected = "")
    .updateButton(rv)
  })

  observeEvent(input$de_file, {
    rv$de_file=input$de_file
    #updateSelectizeInput(session, "example", selected = "")
    .updateButton(rv)
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

  #Take a tour
  observeEvent(input$tour,{
    introjs(session )
  })

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

            res1 = NicheSIG_pipeline (input$species, rv$cond1_file$datapath, input$cutoff, rv$de_file$datapath, input$pctile, invert_DE = FALSE)
            output$results1_active <- DT::renderDataTable({
              .trimResults(res1, active = TRUE)
            },server = TRUE, selection = "single")

            output$results1_inactive <- DT::renderDataTable({
              .trimResults(res1, active = FALSE)
            },server = TRUE, selection = "single")
          })

          output$titleText1 <- renderText(paste("Condition1:",rv$cond1_file$name))

          incProgress(0.4)

          withProgress(message = 'Condition 2',  {
            res2= NicheSIG_pipeline (input$species, rv$cond2_file$datapath, input$cutoff, rv$de_file$datapath, input$pctile, invert_DE = TRUE)
            output$results2_active <- DT::renderDataTable({
              .trimResults(res2, active = TRUE)
            },server = TRUE, selection = "single")

            output$results2_inactive <- DT::renderDataTable({
              .trimResults(res2, active = FALSE)
            },server = TRUE, selection = "single")
          })

          output$titleText2 <- renderText(paste("Condition2:",rv$cond2_file$name))
        })

      }
      , error = function(e){
        showModal(modalDialog(
          title = "Error", e$message,
          footer = modalButton("Dismiss")
        ))
      }
    )

  })

}

#################################
## Running the app (only in interactive environment)
#################################
#if (interactive()) {
  shinyApp(ui, server)
#}
