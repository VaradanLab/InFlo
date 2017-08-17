library(shiny)
library(visNetwork)
library(survival)
require(shinyjs)
require(markdown)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(png)
#library("ReporteRs")
library(shinydashboard)
library(leaflet)

header <- dashboardHeader(title = tags$p("InFlo", style = "font-size: 200%;"),
                          dropdownMenu(type = "messages",
                                       messageItem(
                                         from = "New User",
                                         message = "How do I register?",
                                         icon = icon("question"),
                                         time = "13:45"
                                       ),
                                       messageItem(
                                         from = "Support",
                                         message = "The new server is ready.",
                                         icon = icon("life-ring"),
                                         time = "2014-12-01"
                                       )
                          ),
                          dropdownMenu(type = "notifications",
                                       notificationItem(
                                         text = "5 new users today",
                                         icon("users")
                                       ),
                                       notificationItem(
                                         text = "12 InFlo Running",
                                         icon("truck"),
                                         status = "success"
                                       ),
                                       notificationItem(
                                         text = "Server load at 86%",
                                         icon = icon("exclamation-triangle"),
                                         status = "warning"
                                       )
                          ),
                          dropdownMenu(type = "tasks", badgeStatus = "success",
                                       taskItem(value = 100, color = "green",
                                                "Data Import"
                                       ),
                                       taskItem(value = 100, color = "aqua",
                                                "Data Pre-Processing"
                                       ),
                                       taskItem(value = 75, color = "yellow",
                                                "InFlo"
                                       ),
                                       taskItem(value = 10, color = "red",
                                                "InFlo Post-Processing"
                                       ),
                                       taskItem(value = 80, color = "red",
                                                "Overall project"
                                       )
                          ),
                          titleWidth = 400
)



Sidebar <- dashboardSidebar(
  width = 400,
  fileInput('file5', 'Choose Component File', 
            accept=c('text/tsv',
                     'text/comma-separated-values,text/plain',
                     '.csv')),
  tags$hr(),
  fileInput('file6', 'Choose Interaction File', 
            accept=c('text/tsv',
                     'text/comma-separated-values,text/plain',
                     '.csv')),
  tags$hr(),
  fileInput('file7', 'Choose InFlo File', 
            accept=c('text/tsv',
                     'text/comma-separated-values,text/plain',
                     '.csv')),
  tags$hr(),
  fileInput('file8', 'Choose Clinical File', 
            accept=c('text/tsv',
                     'text/comma-separated-values,text/plain',
                     '.csv')),
  tags$hr(),
  
  div(style="display: inline-block;vertical-align:center;", checkboxInput(inputId="allSamps", label="All Samples in Cohort", value=T)),
  br(),
  tags$hr(),
  div(style="display: inline-block;vertical-align:center;",radioButtons("physics", "Automatic Movement", choices = c("On", "Off"), inline=T)),
  br(),
  tags$hr(),
  div(style="display: inline-block;vertical-align:center;",radioButtons("layout", "Layout", choices = c("Hairball", "Hierarchical"), inline=T)),
  br(),
  tags$hr(),
  actionButton(inputId="neighbors", label="Show Neighbors"),
  actionButton(inputId="remove", label="Remove Selected"),
  br(),
  tags$hr(),
  radioButtons("type", label="Show Nodes That Are", choices=c("Most Deviant", "Most Disparate"), inline=T),
  br(),
  tags$hr()

)


body <- dashboardBody(
  fluidPage(
    navbarPage("Check",id="nav",
               tabPanel("Network",
                        div(class="outer",
                        tags$head(
                              includeCSS("styles.css")
                            ),
                        tags$style(type = "text/css", ".outer {position: fixed; top: 41px; left: 0; right: 0; bottom: 0; overflow: hidden; padding: 0}")
                        ),
                        # Shiny versions prior to 0.11 should use class="modal" instead.
                        mainPanel(width = "100%",height = "auto",
                          fluidPage(
                        tabItem(tabName = "t_item1", class = "active",
                                #box(title = "Network",status = "primary", solidHeader = T, collapsible = T,width ="100%",height = 10,
                          visNetworkOutput("network",width = "100%",height ="1000px")
                          #))
                        ))
               ),
               absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                             draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
                             width = 330, height = "auto",
                             
                             h2("Navigation Panel"),
                             
                             div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("spread1")),
                             br(),
                             div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("Gr")),
                             br(),
                             div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("Gr1")),
                             br(),
                             div(style="display: inline-block;vertical-align:center; width: 250px;",uiOutput("selectSamp")),
                             br(),
                             div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("dis1")),
                             br(),        
                             div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("dis2")),
                             br(),
                             actionButton(inputId="refreshNet", label="Create New Network", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                             actionButton(inputId="revert", label="Original Network"),
                             br(),
                             actionButton(inputId="back", label="<"),
                             actionButton(inputId="forward", label=">"),
                             tags$hr(),
                             actionButton(inputId="plot", label="Refresh Plots"),
                             br()
                             #style = "opacity: 0.92"
               )),
               tabPanel("Data explorer",
                        fluidRow(
                          DT::dataTableOutput('InFlo')
                        )
               ),
               tabPanel("Survival Plot",
                        fluidRow(
                          plotOutput("bar", width="50%"),
                          br(),
                          plotOutput('allSurvivalCurve'),
                          br(),
                          plotOutput('G1survivalCurve'),
                          br(),
                          plotOutput('G2survivalCurve')
                        )
               )
               
               
               
               
      )
  )
)
  
  
#     mainPanel(
#                div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("Gr")),
#                div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("Gr1")),
#                #br(),
#                div(style="display: inline-block;vertical-align:center; width: 250px;",uiOutput("selectSamp")),
#                #div(style="display: inline-block;vertical-align:center;", checkboxInput(inputId="allSamps", label="All Samples in Cohort", value=T)),
#                #br(),
#                #div(style="display: inline-block;vertical-align:center;",radioButtons("physics", "Automatic Movement", choices = c("On", "Off"), inline=T)),
#                div(style="display: inline-block;vertical-align:top; width: 100px;",HTML("<br>")),
#                #div(style="display: inline-block;vertical-align:center;",radioButtons("layout", "Layout", choices = c("Hairball", "Hierarchical"), inline=T)),
#                #br(),
#                visNetworkOutput("network", height='600px'),
#                actionButton(inputId="neighbors", label="Show Neighbors"),
#                actionButton(inputId="remove", label="Remove Selected"),
#                #br(),
#                #br(),
#                #br(),
#                #radioButtons("type", label="Show Nodes That Are", choices=c("Most Deviant", "Most Disparate"), inline=T),
#                # div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("spread1")),
#                # br(),
#                # div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("dis1")),
#                # div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("dis2")),
#                # br(),
#                # br(),
#                # actionButton(inputId="refreshNet", label="Create New Network", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
#                # actionButton(inputId="revert", label="Original Network"),
#                #
#                # br(),
#                # br(),
#                # actionButton(inputId="back", label="<"),
#                # actionButton(inputId="forward", label=">"),
#                #
#                # tags$hr(),
#                # actionButton(inputId="plot", label="Refresh Plots"),
#                # br(),
#                plotOutput("bar", width="50%"),
# 
#                plotOutput('allSurvivalCurve'),
#                plotOutput('G1survivalCurve'),
#                plotOutput('G2survivalCurve'),
#                br(),
#                br(),
#                tags$hr(),
#                uiOutput("selectNames"),
#                br(),
#                actionButton(inputId="refreshInflo", label="Refresh Table"),
#                br(),
#                dataTableOutput('InFlo'),
#                tags$hr()
#                )
#       )
# )
    
dashboardPage(
  header,
  Sidebar,
  #dashboardSidebar(disable = TRUE),
  body
)
#   # tags$head(tags$style(HTML('
#   #                           /* main sidebar */
#   #                           .skin-blue .main-sidebar {
#   #                           background-color: #c9faff;
#   #                           }'
# #)))
#   )


# shinyUI(
#   fluidPage(
#     
#     ############################################################ Header #####################################################
#     titlePanel(
#       img(src="img1.png", height = "100", width = "100", align = 'left')
#     ),
#     titlePanel(h1("InFlo", align = "center",style="color:blue")),
#     titlePanel(h3("An integrative systems biology framework for characterizing activities of complex signaling networks in individual patient samples using multi-omics profiling", 
#                   align = "center",style="color:blue")
#     ),
#     br(),
#     tags$hr(),
#     ################################################################################################################
#     ########################################Side Absolute Panel #######################################################
#     absolutePanel(
#       bottom = 20, right = 20, width = 300,
#       draggable = TRUE,
#       fixed = F,
#         div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("spread1")),
#         br(),
#         div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("dis1")),
#         div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("dis2")),
#         br(),
#         br(),
#         actionButton(inputId="refreshNet", label="Create New Network", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
#         actionButton(inputId="revert", label="Original Network"),
#         
#         br(),
#         br(),
#         actionButton(inputId="back", label="<"),
#         actionButton(inputId="forward", label=">"),
#         
#         tags$hr(),
#         actionButton(inputId="plot", label="Refresh Plots"),
#         br(),
#       style = "opacity: 0.92"
#         ),
#     ####################################################################################################################
#     ########################################################Side Panel##############################################
#     sidebarLayout(
#       sidebarPanel(
#         #fileInput('file1', 'Choose Sample File',
#         #          accept=c('text/tsv',
#         #                   'text/comma-separated-values,text/plain',
#         #                   '.csv')),
#         #tags$hr(),
#         #fileInput('file2', 'Choose Expression File',
#         #          accept=c('text/tsv',
#         #                   'text/comma-separated-values,text/plain',
#         #                   '.csv')),
#         #tags$hr(),
#         #fileInput('file3', 'Choose CopyNumber File',
#         #          accept=c('text/tsv',
#         #                   'text/comma-separated-values,text/plain',
#         #                   '.csv')),
#         #tags$hr(),
#         #fileInput('file4', 'Choose Pathway File',
#         #          accept=c('text/tsv',
#         #                   'text/comma-separated-values,text/plain',
#         #                   '.csv')),
#         #tags$hr(),
#         
#         fileInput('file5', 'Choose Component File', 
#                   accept=c('text/tsv',
#                            'text/comma-separated-values,text/plain',
#                            '.csv')),
#         tags$hr(),
#         fileInput('file6', 'Choose Interaction File', 
#                  accept=c('text/tsv',
#                           'text/comma-separated-values,text/plain',
#                           '.csv')),
#         tags$hr(),
#         fileInput('file7', 'Choose InFlo File', 
#                   accept=c('text/tsv',
#                            'text/comma-separated-values,text/plain',
#                            '.csv')),
#         tags$hr(),
#         fileInput('file8', 'Choose Clinical File', 
#                   accept=c('text/tsv',
#                            'text/comma-separated-values,text/plain',
#                            '.csv')),
#         tags$hr(),
#         
#         div(style="display: inline-block;vertical-align:center;", checkboxInput(inputId="allSamps", label="All Samples in Cohort", value=T)),
#         br(),
#         tags$hr(),
#         div(style="display: inline-block;vertical-align:center;",radioButtons("physics", "Automatic Movement", choices = c("On", "Off"), inline=T)),
#         br(),
#         tags$hr(),
#         div(style="display: inline-block;vertical-align:center;",radioButtons("layout", "Layout", choices = c("Hairball", "Hierarchical"), inline=T)),
#         br(),
#         tags$hr(),
#         actionButton(inputId="neighbors", label="Show Neighbors"),
#         actionButton(inputId="remove", label="Remove Selected"),
#         br(),
#         tags$hr(),
#         radioButtons("type", label="Show Nodes That Are", choices=c("Most Deviant", "Most Disparate"), inline=T),
#         br(),
#         tags$hr(),
#         #radioButtons('sep', 'Expression DataType',
#         #             c(RNASeqV2 = T,
#         #               Microarray = T),
#         #             'F'),
#         # radioButtons('quote', 'Quote',
#         #              c(None='',
#         #                'Double Quote'='"',
#         #                'Single Quote'="'"),
#         #              '"'),
#         #tags$hr(),
#         #actionButton("goButton", "Go!"),
#         tags$hr()
#         # checkboxGroupInput('show_vars', 'Columns in Results to show:',
#         #                    names(header), selected = names(header))
#         #actionButton("goButton","Go!"),
#       ),
#       ####################################################################################################################
#       #################################### Main Panel #################################################################
#       mainPanel(
#          div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("Gr")),
#          div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("Gr1")),
#          #br(),
#          div(style="display: inline-block;vertical-align:center; width: 250px;",uiOutput("selectSamp")), 
#          #div(style="display: inline-block;vertical-align:center;", checkboxInput(inputId="allSamps", label="All Samples in Cohort", value=T)),
#          #br(),
#          #div(style="display: inline-block;vertical-align:center;",radioButtons("physics", "Automatic Movement", choices = c("On", "Off"), inline=T)),
#          div(style="display: inline-block;vertical-align:top; width: 100px;",HTML("<br>")),
#          #div(style="display: inline-block;vertical-align:center;",radioButtons("layout", "Layout", choices = c("Hairball", "Hierarchical"), inline=T)),
#          #br(),
#          visNetworkOutput("network", height='600px'),
#          actionButton(inputId="neighbors", label="Show Neighbors"),
#          actionButton(inputId="remove", label="Remove Selected"),
#          #br(),
#          #br(),
#          #br(),
#          #radioButtons("type", label="Show Nodes That Are", choices=c("Most Deviant", "Most Disparate"), inline=T),
#          # div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("spread1")),
#          # br(),
#          # div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("dis1")),
#          # div(style="display: inline-block;vertical-align:center; width: 250px;", uiOutput("dis2")),
#          # br(),
#          # br(),
#          # actionButton(inputId="refreshNet", label="Create New Network", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
#          # actionButton(inputId="revert", label="Original Network"),
#          # 
#          # br(),
#          # br(),
#          # actionButton(inputId="back", label="<"),
#          # actionButton(inputId="forward", label=">"),
#          # 
#          # tags$hr(),
#          # actionButton(inputId="plot", label="Refresh Plots"),
#          # br(),
#          plotOutput("bar", width="50%"),
#          
#          plotOutput('allSurvivalCurve'),
#          plotOutput('G1survivalCurve'),
#          plotOutput('G2survivalCurve'),
#          br(),
#          br(),
#          tags$hr(),
#          uiOutput("selectNames"),
#          br(),
#          actionButton(inputId="refreshInflo", label="Refresh Table"),
#          br(),
#          dataTableOutput('InFlo'),
#          tags$hr()
#          
#   )
# )
# )
# )
