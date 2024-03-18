# Load library
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
      # selectInput("GeneSplicer", "GeneSplicer cutoff:", choices = c("50", "100", "200")),

      # Choose SQUIRLS cutoff
      sliderInput("SQUIRLS", "SQUIRLS cutoff:", min = 0, max = 1, value = 0.9, step = 0.1),

      # Choose MMSplice cutoff
      sliderInput("MMSplice", "MMSplice cutoff:", min = 0, max = 2, value = 0.5, step = 0.5),

      # Choose Pangolin cutoff
      #sliderInput("Pangolin", "Pangolin cutoff:", min = 0, max = 1, value = 0.2, step = 0.1),
      
      # Choose what to plot
      selectInput("plot", "FDR plot:", choices = c("SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", "SpliceAI_DS_DL", "MaxEntScan_alt", "SQUIRLS", "mmsplice_delta_logit_psi"))
      
  ),

  # Create a main panel to display table of results
  mainPanel(

    # Table of splice variants and prediction scores
    box(style='width:1000px;overflow-x: scroll; overflow-y: scroll;',
    DT::dataTableOutput("splice_table")
    ),
    
    # Table of performance metrics
    DT::dataTableOutput("metrics_table"),
     
    # FP-FN trade-off plot 
    plotOutput(outputId = "main_plot", height = "300px")
  )
)

# Define server
server <- function(input, output) {

  # Filter data based on selections
  output$splice_table <- DT::renderDataTable(DT::datatable({
    #file <- paste("/data/jess_tmp/fh/Rshiny/splice/", input$worklist, ".csv", sep="")
    file <- paste(input$worklist, ".csv", sep="")
    
    data <- read.csv(file)
    
    if (input$SpliceAI != "All") {
      data <- subset(data, SpliceAI_DS_AG > input$SpliceAI | SpliceAI_DS_AL > input$SpliceAI | SpliceAI_DS_DG > input$SpliceAI | SpliceAI_DS_DL > input$SpliceAI)
    }
    
    if (input$MES != "All") {
      data <- subset(data, (MaxEntScan_diff < 0 & MaxEntScan_alt > input$MES) | (MaxEntScan_diff > 0 & MaxEntScan_alt < input$MES))
    }
    
    # if (input$GeneSplicer != "All") {
    #   data <- subset(data, GeneSplicer_score > input$GeneSplicer)
    # }
    
    #if (input$MMSplice != "All") {
      #data <- subset(data, (mmsplice_delta_logit_psi > input$MMSplice) | (mmsplice_delta_logit_psi < (input$MMSplice * -1)))
    #}
    
    if (input$MMSplice != "All") {
      data <- rbind(data[as.numeric(data$mmsplice_delta_logit_psi) < (-1 * input$MMSplice)], data[as.numeric(data$mmsplice_delta_logit_psi) > input$MMSplice] )
    }
     
    #if (input$Pangolin != "All") {
      #data <- subset(data, Pangolin_score_change_1 > input$Pangolin | Pangolin_score_change_2 > input$Pangolin)
    #}
    
    if (input$SQUIRLS != "All") {
      data <- subset(data, SQUIRLS > input$SQUIRLS)
    }
    
    data
  }))

   # Performance metrics TP, FP, FDR = FP / (FP + TP)
   output$metrics_table <- DT::renderDataTable(DT::datatable({
    
    # Get metrics
    TP <- nrow(data[data$Type == "TP",])
    FP <- nrow(data[data$Type == "FP",])
    FDR <- FP / (FP + TP)
         
    # create matrix with 4 columns and 4 rows
    table= matrix(c(TP, FP, FDR), ncol=3, byrow=TRUE)
      
    # specify the column names and row names of matrix
    colnames(table) <- c('TP','FP','FDR')
 
    # assign to table
    metrics=as.table(table)
     
    metrics
  }))
    
    # FP-FN trade-off plot https://stackoverflow.com/questions/70841834/false-positive-vs-false-negative-trade-off-plot    https://stackoverflow.com/questions/6939136/how-to-overlay-density-plots-in-r 
    output$main_plot <- renderPlot({
        # Get data
        column <- input$plot
        TP_dat <- data[data$Type == 'TP', ]
        FP_dat <- data[data$Type == 'FP', ]
         
        # Plot
        plot(density(TP_dat$column))
        lines(density(FP_dat$column))
  })

}

# Run app
shinyApp(ui = ui, server = server)
