# Define UI for application that plots random distributions 
if (!require("shinydashboard")){
  install.packages("shinydashboard")
  library(shinydashboard)
}else{
  library(shinydashboard)
}

shinyUI(
  dashboardPage(
    dashboardHeader(
      title = "",
      titleWidth = "25%"
    ),
    dashboardSidebar(
      width = "25%",
      sidebarUserPanel("INRA TOULOUSE 2018",
                       subtitle = "Please use the «Add file» button below to upload a file",
                       image = "images/inra.png"
      ),
      fluidPage(
      actionButton("addSelect", "Add file", icon = icon("plus", class="fas fa-upload", 
                                                        lib = "font-awesome")),
      tags$hr(),
      span(id = 'placeholder'),
      tags$style("#addSelect{cursor: pointer; opacity: 0.8;}")
      )
    ),
    # body
    dashboardBody(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "www/ressources/astyles.css")
      ),
      tabsetPanel(
        tabPanel("Venn Diagram", tags$hr(), plotOutput("outPlot", width = "auto", height = "850px")),
        tabPanel("Upset Plot", tags$hr(), plotOutput("upsetPlot", width = "auto", height = "850px")),
        tabPanel("Summary", tags$hr(), verbatimTextOutput("summary"),
                 tags$head(tags$style("#summary{ overflow-y: scroll; max-height: 850px; }"))
                 ),
        tabPanel("ECs Tab", tags$hr(), 
           fluidRow(
             column(2,
             actionButton("addSelectTab", "Select a file", icon=icon("plus", class="fas fa-upload", 
                                                                       lib = "font-awesome")),
             tags$style("#addSelectTab{cursor: pointer;}")
             ),
             column(6),
             column(4, 
                    textInput(inputId = "inputECnumber", label = "Enter an EC Number", value = "")
                  ),
             tags$style(type='text/css', "#addSelectTab {width:60%; margin-left: 25px; top:0px}"),
             tags$style(type='text/css', "#inputECnumber {width:60%; margin-right: 25px; top:0px}")
           ),
         fluidRow(
           span(id = 'placeholderTab')
         ),
         fluidPage(
          tableOutput("ectab"),
         tags$head(tags$style("#ectab{ overflow-y: scroll; max-height: 650px; }"))
         )
        ),
        tabPanel("Heatmap", tags$hr(), 
                 fluidRow(
                   column(2,
                          actionButton("addSelectHeat", "Select a file", icon=icon("plus", class="fas fa-upload", 
                                                                                  lib = "font-awesome")),
                          tags$style("#addSelectHeat{cursor: pointer;}")
                   ),
                   tags$style(type='text/css', "#addSelectHeat {width:60%; margin-left: 25px; top:0px}")
                 ),
                 fluidRow(
                   span(id = 'placeholderHeat')
                 ),
                 plotOutput("heatmapPlot", width = "auto", height = "850px")
                 # imageOutput("heatmapPlot", width = "auto", height = "850px")
        )
      )
    )
  )
)