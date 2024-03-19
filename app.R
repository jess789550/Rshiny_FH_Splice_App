# Load libraries
library(shiny)
library(DT)
library(shinydashboard) # for box()
library(ggplot2)

# Define UI: https://shiny.posit.co/r/gallery/widgets/basic-datatable/
ui <- fluidPage(
  # Title
  titlePanel("FH Splice Site Prediction Results"),
  
  # Create a sidebar panel in the UI for adjusting parameters
  sidebarPanel(
    
    # Allow user to select worklist
    selectInput("worklist", "Choose a worklist:",
                choices = c("2004442", "2005265", "2005267", "2005745", "2324533",
                            "2400720", "2322015", "2127291", "2226732", "2327211",
                            "2330804", "2331473", "mix")),
    
    # Choose SpliceAI cutoff
    sliderInput("SpliceAI", "SpliceAI cutoff:", min = 0, max = 1, value = 0.3, step = 0.1),
    
    # Choose MES cutoff
    sliderInput("MES", "MES cutoff:", min = 0, max = 10, value = 6.2, step = 0.1),
    
    # Choose GeneSplicer cutoff
    # selectInput("GeneSplicer", "GeneSplicer cutoff:", choices = c("50", "100", "200", "None")),
    
    # Choose SQUIRLS cutoff
    sliderInput("SQUIRLS", "SQUIRLS cutoff:", min = 0, max = 1, value = 0.9, step = 0.1),
    
    # Choose MMSplice cutoff
    sliderInput("MMSplice", "MMSplice cutoff:", min = 0, max = 2, value = 0.5, step = 0.5),
    
    # Choose Pangolin cutoff
    #sliderInput("Pangolin", "Pangolin cutoff:", min = 0, max = 1, value = 0.2, step = 0.1),
    
    # Choose what to plot
    selectInput("plot", "FDR plot:", choices = c("SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", 
                                                 "SpliceAI_DS_DL", "MaxEntScan_alt", "SQUIRLS", "mmsplice_delta_logit_psi")),
    
    # Add button for user to press to initiate run
    actionButton(inputId = "Submit", label = "Submit")
  ),
  
  # Create a main panel to display table of results
  mainPanel(
    
    # Add description from README.md
    uiOutput('Description'),
    br(),
    
    # Table of performance metrics
    h2("Table of performance metrics"),
    DT::dataTableOutput("metrics_table"),
    br(),
    
    # TP-FP trade-off plot 
    h2("TP-FP trade-off plot"),
    plotOutput(outputId = "main_plot", height = "300px"),
    
    # Table of splice variants and prediction scores
    h2("Table of splice variants and prediction scores"),
    box(style='width:1000px;overflow-x: scroll; overflow-y: scroll;',
        DT::dataTableOutput("splice_table")
    ),
    br()
  )
)

# Define server
server <- function(input, output) {
  
  # Get description from README.md
  output$Description <- renderUI({
    rawText <- readLines('README.md') # get raw text
    
    # split the text into a list of character vectors
    #   Each element in the list contains one line
    splitText <- stringi::stri_split(str = rawText, regex = '\\n')
    
    # wrap a paragraph tag around each element in the list
    replacedText <- lapply(splitText, p)
    
    return(replacedText)
  })
  
  observeEvent(input$Submit, {
      
      # Filter data based on selections
      file <- paste(input$worklist, ".csv", sep="")
      
      data <- read.csv(file)
      
      if (input$SpliceAI != 0) {
        data <- subset(data, SpliceAI_DS_AG > input$SpliceAI | SpliceAI_DS_AL > input$SpliceAI | 
                         SpliceAI_DS_DG > input$SpliceAI | SpliceAI_DS_DL > input$SpliceAI)
      }
      
      if (input$MES != 0) {
        dataFiltered<-data[which(data$MaxEntScan_diff!='-'),]
        data <- subset(dataFiltered, (MaxEntScan_diff < 0 & MaxEntScan_alt > input$MES) | 
                         (MaxEntScan_diff > 0 & MaxEntScan_alt < input$MES))
      }
      
      # if (input$GeneSplicer != "None") {
      #   data <- subset(data, GeneSplicer_score > input$GeneSplicer)
      # }
      
      #if (input$MMSplice != 0) {
      #data <- subset(data, (mmsplice_delta_logit_psi > input$MMSplice) | (mmsplice_delta_logit_psi < (input$MMSplice * -1)))
      #}
      
      if (input$MMSplice != 0) {
        dataFiltered<-data[which(data$mmsplice_delta_logit_psi!='-'),]
        data <- rbind(data[as.numeric(dataFiltered$mmsplice_delta_logit_psi) < (-1 * input$MMSplice),], 
                      data[as.numeric(dataFiltered$mmsplice_delta_logit_psi) > input$MMSplice,])
      }
      
      #if (input$Pangolin != 0) {
      #data <- subset(data, Pangolin_score_change_1 > input$Pangolin | Pangolin_score_change_2 > input$Pangolin)
      #}
      
      if (input$SQUIRLS != 0) {
        data <- subset(data, SQUIRLS > input$SQUIRLS)
      }
      dataDebug<<-data
    
      # Show table of filtered data
      output$splice_table <- DT::renderDataTable(DT::datatable({
        data
      }))
      
      # Performance metrics TP, FP, FDR = FP / (FP + TP)
      output$metrics_table <- DT::renderDataTable(DT::datatable({
        
        # Get metrics
        TP <- nrow(data[data$Type == "TP",])
        FP <- nrow(data[data$Type == "FP",])
        FDR <- FP / (FP + TP)
        
        # create matrix 
        table= matrix(c(TP, FP, FDR), ncol=3, byrow=TRUE)
        
        # specify the column names and row names of matrix
        colnames(table) <- c('TP','FP','FDR')
        
        # assign to table
        metrics=as.data.frame(table)
        
        metrics
        metricsDebug<<-metrics
      }))
      
      # Get data
      #TP_dat <- data[data$Type == 'TP', ]
      #FP_dat <- data[data$Type == 'FP', ]
      
      # Plot input
      plot_func <- function(column) {
        renderPlot({
          # plot(density(TP_dat[,column][!is.na(TP_dat[,column])]))
          # lines(density(FP_dat[,column][!is.na(FP_dat[,column])]))
          column <- sym(column)
          if (column == "MaxEntScan_alt") {
            dataFiltered<-data[which(data$MaxEntScan_diff!='-'),]
            ggplot(dataFiltered, aes(x = !!column, fill = Type)) + geom_density(alpha = 0.5)
          } else if (column =="mmsplice_delta_logit_psi") {
            dataFiltered<-data[which(data$mmsplice_delta_logit_psi!='-'),]
            ggplot(dataFiltered, aes(x = !!column, fill = Type)) + geom_density(alpha = 0.5)
          } else{
            ggplot(data, aes(x = !!column, fill = Type)) + geom_density(alpha = 0.5)
          }
        })
      }
        #column <- input$plot
      
      # TP-FP trade-off plot 
      # https://stackoverflow.com/questions/70841834/false-positive-vs-false-negative-trade-off-plot   
      # https://stackoverflow.com/questions/6939136/how-to-overlay-density-plots-in-r 
      output$main_plot <- tryCatch(
        {
          plot_func(input$plot)
        },
        error = function(cond) {
          message("Sorry the TP-FP trade-off plot cannot be produced.")
          message("This could be due to no TP/FP in your selection.")
          message("Here's the original error message:")
          message(conditionMessage(cond))
          # Choose a return value in case of error
          NA
        },
        warning = function(cond) {
          message("Sorry the TP-FP trade-off plot cannot be produced.")
          message("This could be due to no TP/FP in your selection.")
          message("Here's the original warning message:")
          message(conditionMessage(cond))
          # Choose a return value in case of warning
          NULL
        }
      )
    })
  }
  
  # Run app
  shinyApp(ui = ui, server = server)
  
