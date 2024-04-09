# DEPLOYMENT VERSION #

##### Load libraries ######
library(shiny)
library(DT)
library(shinydashboard) # for box()
library(ggplot2)

##### Define UI: https://shiny.posit.co/r/gallery/widgets/basic-datatable/ #####
ui <- dashboardPage(
  ### Title ##
  dashboardHeader(
    title="FH Splice Site Prediction Results"
  ),
  
  ### Create a sidebar panel in the UI for adjusting parameters ###
  dashboardSidebar(
    
    # Allow user to select worklist
    selectInput("worklist", "Choose a worklist:",
                choices = c("2004442", "2005265", "2005267", "2005745", "2324533",
                            "2400720", "2322015", "2127291", "2226732", "2327211",
                            "2330804", "2331473", "mix")),
    
    # Choose SpliceAI cutoff
    sliderInput("SpliceAI", "SpliceAI cutoff:", min = 0, max = 1, value = 0.3, step = 0.1),
    
    # Choose MES cutoff
    # sliderInput("MES", "MES cutoff:", min = 0, max = 10, value = 6.2, step = 0.1),
    selectInput("MES", "MES cutoff:", choices = c("Low", "None", "High")),
    
    # Choose GeneSplicer cutoff
    # selectInput("GeneSplicer", "GeneSplicer cutoff:", choices = c("50", "100", "200", "None")),
    
    # Choose SQUIRLS cutoff
    sliderInput("SQUIRLS", "SQUIRLS cutoff:", min = 0, max = 1, value = 0.5, step = 0.1),
    
    # Choose MMSplice cutoff
    sliderInput("MMSplice", "MMSplice cutoff:", min = 0, max = 2, value = 0.5, step = 0.5),
    
    # Choose Pangolin cutoff
    #sliderInput("Pangolin", "Pangolin cutoff:", min = 0, max = 1, value = 0.2, step = 0.1),
    
    # Choose gnomAD cutoff
    # sliderInput("gnomAD", "gnomAD allele frequency:", min = 0, max = 1, value = 0, step = 0.1),
    selectInput("gnomAD", "gnomAD allele frequency cutoff:", choices = c("None", "0", "0.0001", "0.0002", 
                                                                         "0.002","0.005")),
    
    # Choose what to plot
    selectInput("plot", "FDR plot:", choices = c("SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", 
                                                 "SpliceAI_DS_DL", "MaxEntScan_alt", "SQUIRLS", "mmsplice_delta_logit_psi")),
    
    # Add button for user to press to initiate run
    actionButton(inputId = "Submit", label = "Submit")
  ),
  
  ### Create a main panel to display table of results ###
  dashboardBody(
    
    tabsetPanel(type = "tabs",
                tabPanel("Description" ,
                         # Add description from README.md
                         uiOutput('Description'),
                ),
                tabPanel("Splice variants and prediction scores",
                         # Table of splice variants and prediction scores
                         h2("Table of splice variants and prediction scores"),
                         fluidRow(box(style='width:1000px;',
                                      fluidRow(
                                        column(2, selectizeInput(inputId="col_file_id", label="File:", choices="All")),
                                        column(2, selectizeInput(inputId="col_CHROM", label="Chromosome:", choices="All")),
                                        column(2,selectizeInput(inputId="col_POS", label="Position:", choices="All")),
                                        column(2,selectizeInput(inputId="col_REF", label="Reference allele:", choices="All")),
                                        column(2,selectizeInput(inputId="col_ALT", label="Alternative allele:", choices="All")),
                                        column(2,selectizeInput(inputId="col_SYMBOL", label="Gene symbol:", choices="All")),
                                      ),
                                      fluidRow(
                                        column(2, selectizeInput(inputId="col_HGVSc", label="HGVSc nomenclature:", choices="All")),
                                        column(2, selectizeInput(inputId="col_gnomAD_AF", label="gnomAD allele frequency:", choices="All")),
                                        column(2, selectizeInput(inputId="col_SpliceAI_DS_AG", label="SpliceAI Acceptor Gain:", choices="All")),
                                        column(2, selectizeInput(inputId="col_SpliceAI_DS_AL", label="SpliceAI Acceptor Loss:", choices="All")),
                                        column(2, selectizeInput(inputId="col_SpliceAI_DS_DG", label="SpliceAI Donor Gain:", choices="All")),
                                        column(2, selectizeInput(inputId="col_SpliceAI_DS_DL", label="SpliceAI Donor Loss:", choices="All")),
                                      ),
                                      fluidRow(
                                        column(2, selectizeInput(inputId="col_mmsplice_delta_logit_psi", label="MMSplice:", choices="All")),
                                        column(2, selectizeInput(inputId="col_MaxEntScan_alt", label="MES alt:", choices="All")),
                                        column(2, selectizeInput(inputId="col_MaxEntScan_diff", label="MES diff:", choices="All")),
                                        column(2, selectizeInput(inputId="col_SQUIRLS", label="SQUIRLS:", choices="All")),
                                        column(2, selectizeInput(inputId="col_Type", label="True positive or false positive:", choices="All")),
                                        column(2, selectizeInput(inputId="col_QUAL", label="VCF Quality score:", choices="All")),
                                      ),
                                      fluidRow(
                                        column(2, actionButton(inputId = "Filter", label = "Filter"))
                                      ))),
                         br(),
                         fluidRow(box(style='width:800px;overflow-x: scroll; overflow-y: scroll;',
                                      DT::dataTableOutput("splice_table"),))
                ),
                tabPanel("Performance metrics",
                         # Table of performance metrics
                         h2("Table of performance metrics"),
                         DT::dataTableOutput("metrics_table")
                ),
                tabPanel("TP-FP trade-off plot",
                         # TP-FP trade-off plot 
                         h2("TP-FP trade-off plot for all data"),
                         plotOutput(outputId = "main_plot", height = "300px"),
                         br(),
                         h2("TP-FP trade-off plot for filtered data"),
                         plotOutput(outputId = "filtered_plot", height = "300px")
                )
    )
  )
)

##### Define server #####
server <- function(input, output, session) {
  
  ### Get description from README.md ###
  output$Description <- renderUI({
    rawText <- readLines('README.md') # get raw text
    
    # split the text into a list of character vectors
    #   Each element in the list contains one line
    splitText <- stringi::stri_split(str = rawText, regex = '\\n')
    
    # wrap a paragraph tag around each element in the list
    replacedText <- lapply(splitText, p)
    
    return(replacedText)
  })
  
  ### On click of submit button produce tables and plot ###
  observeEvent(input$Submit, {
    
    # Read worklist splice site prediction results
    file <- paste(input$worklist, ".csv", sep="")
    original_data <- read.csv(file)  # need original_data downstream
    data <- original_data  # for filtering
    dataDebug_Intro<<-data
    # Update column filters 
    # https://stackoverflow.com/questions/46346917/update-shinys-selectinput-dropdown-with-new-values-after-uploading-new-data-u
    updateSelectizeInput(session, "col_file_id",
                         choices = c("All", unique(as.character(data$file_id))))
    
    updateSelectizeInput(session, "col_CHROM",
                         choices = c("All", unique(as.character(data$CHROM))))
    
    # updateSelectizeInput(session, "col_POS",
    #                   choices = c("All", unique(as.character(data$POS))))
    
    updateSelectizeInput(session, "col_REF",
                         choices = c("All", unique(as.character(data$REF))))
    
    updateSelectizeInput(session, "col_ALT",
                         choices = c("All", unique(as.character(data$ALT))))
    
    # updateSelectizeInput(session, 'col_SYMBOL', 
    #                      choices = c("All", unique(as.character(data$SYMBOL))))
    # 
    # updateSelectizeInput(session, "col_HGVSc",
    #                   choices = c("All", unique(as.character(data$HGVSc))))
    
    updateSelectizeInput(session, "col_gnomAD_AF",
                         choices = c("All", unique(as.character(data$gnomAD_AF))))
    
    updateSelectizeInput(session, "col_SpliceAI_DS_AG",
                         choices = c("All", unique(as.character(data$SpliceAI_DS_AG))))
    
    updateSelectizeInput(session, "col_SpliceAI_DS_AL",
                         choices = c("All", unique(as.character(data$SpliceAI_DS_AL))))
    
    updateSelectizeInput(session, "col_SpliceAI_DS_DG",
                         choices = c("All", unique(as.character(data$SpliceAI_DS_DG))))
    
    updateSelectizeInput(session, "col_SpliceAI_DS_DL",
                         choices = c("All", unique(as.character(data$SpliceAI_DS_DL))))
    
    # updateSelectizeInput(session, "col_mmsplice_delta_logit_psi",
    #                   choices = c("All", unique(as.character(data$mmsplice_delta_logit_psi))))
    
    updateSelectizeInput(session, "col_MaxEntScan_alt",
                         choices = c("All", unique(as.character(data$MaxEntScan_alt))))
    
    updateSelectizeInput(session, "col_MaxEntScan_diff",
                         choices = c("All", unique(as.character(data$MaxEntScan_diff))))
    
    # updateSelectizeInput(session, "col_SQUIRLS",
    #                   choices = c("All", unique(as.character(data$SQUIRLS))))
    
    updateSelectizeInput(session, "col_Type",
                         choices = c("All", unique(as.character(data$Type))))
    
    # updateSelectizeInput(session, "col_QUAL",
    #                   choices = c("All", unique(as.character(data$QUAL))))
    
    # Filter SpliceAI results
    if (input$SpliceAI != 0) {
      data <- subset(data, SpliceAI_DS_AG > input$SpliceAI | SpliceAI_DS_AL > input$SpliceAI | 
                       SpliceAI_DS_DG > input$SpliceAI | SpliceAI_DS_DL > input$SpliceAI)
    }
    
    # Filter MES results
    dataFilteredMES<-data[which(data$MaxEntScan_diff!='-'),]
    dataFilteredMES$MaxEntScan_alt<-as.numeric(dataFilteredMES$MaxEntScan_alt)
    
    if (input$MES == "Low") {
      data <- subset(dataFilteredMES, (MaxEntScan_diff < 0 & MaxEntScan_alt > 6.2) | 
                       (MaxEntScan_diff > 0 & MaxEntScan_alt < 8.5))
    } else if (input$MES == "High") {
      data <- subset(dataFilteredMES, (MaxEntScan_diff < 0 & MaxEntScan_alt > 8.5) | 
                       (MaxEntScan_diff > 0 & MaxEntScan_alt < 6.2))
    }
    
    # if (input$GeneSplicer != "None") {
    #   data <- subset(data, GeneSplicer_score > input$GeneSplicer)
    # }
    
    #if (input$MMSplice != 0) {
    #data <- subset(data, (mmsplice_delta_logit_psi > input$MMSplice) | (mmsplice_delta_logit_psi < (input$MMSplice * -1)))
    #}
    
    # Filter MMSplice results
    if (input$MMSplice != 0) {
      dataFilteredMMSplice <- data[which(data$mmsplice_delta_logit_psi!='-'),]
      dataFilteredMMSplice$mmsplice_delta_logit_psi <- as.numeric(dataFilteredMMSplice$mmsplice_delta_logit_psi)
      data <- rbind(data[as.numeric(dataFilteredMMSplice$mmsplice_delta_logit_psi) < (-1 * input$MMSplice),], 
                    data[as.numeric(dataFilteredMMSplice$mmsplice_delta_logit_psi) > input$MMSplice,])
    }
    
    #if (input$Pangolin != 0) {
    #data <- subset(data, Pangolin_score_change_1 > input$Pangolin | Pangolin_score_change_2 > input$Pangolin)
    #}
    
    # Filter SQUIRLS results
    if (input$SQUIRLS != 0) {
      data <- subset(data, SQUIRLS > input$SQUIRLS)
    }
    
    # Filter gnomAD allele frequency
    if (input$gnomAD != "None") {
      freq <- as.numeric(input$gnomAD)
      dataFilteredgnomAD <- data[which(data$gnomAD_AF!='-'),]
      dataFilteredgnomAD$gnomAD_AF <- as.numeric(dataFilteredgnomAD$gnomAD_AF)
      data <- subset(dataFilteredgnomAD, gnomAD_AF <= freq)
    }
    
    filtered_data <<- data
    
    ### Show table of filtered data ###
    output$splice_table <- DT::renderDataTable(DT::datatable({
      data
    }))
    
    ### Performance metrics TP, FP, FDR = FP / (FP + TP) ###
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
    
    # Plot input function
    # column <- input$plot
    plot_func <- function(dataset, column) {
      renderPlot({
        # plot(density(TP_dat[,column][!is.na(TP_dat[,column])]))
        # lines(density(FP_dat[,column][!is.na(FP_dat[,column])]))
        column <- sym(column)
        if (column == "MaxEntScan_alt") {
          dataFilteredMES<-dataset[which(dataset$MaxEntScan_diff!='-'),]
          dataFilteredMES$MaxEntScan_alt<-as.numeric(dataFilteredMES$MaxEntScan_alt)
          ggplot(dataFilteredMES, aes(x = !!column, fill = Type)) + geom_density(alpha = 0.5)
        } else if (column =="mmsplice_delta_logit_psi") {
          dataFilteredMMSplice<-dataset[which(dataset$mmsplice_delta_logit_psi!='-'),]
          dataFilteredMMSplice$mmsplice_delta_logit_psi<-as.numeric(dataFilteredMMSplice$mmsplice_delta_logit_psi)
          ggplot(dataFilteredMMSplice, aes(x = !!column, fill = Type)) + geom_density(alpha = 0.5)
        } else{
          ggplot(dataset, aes(x = !!column, fill = Type)) + geom_density(alpha = 0.5)
        }
      })
    }
    
    ### TP-FP trade-off plot ###
    output$main_plot <- tryCatch(
      {
        plot_func(original_data, input$plot)
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
    
    ### TP-FP trade-off plot filtered ### 
    # https://stackoverflow.com/questions/70841834/false-positive-vs-false-negative-trade-off-plot   
    # https://stackoverflow.com/questions/6939136/how-to-overlay-density-plots-in-r 
    output$filtered_plot <- tryCatch(
      {
        plot_func(data, input$plot)
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
  
  ### Filter table of splice variants and prediction scores ###
  observeEvent(input$Filter, {
    # Filter by column values
    if (input$col_file_id != "All") {
      filtered_data <- filtered_data[filtered_data$file_id == input$col_file_id,]
    }
    
    if (input$col_CHROM!= "All") {
      filtered_data <- filtered_data[filtered_data$CHROM == input$col_CHROM,]
    }
    
    if (input$col_POS!= "All") {
      filtered_data <- filtered_data[filtered_data$POS == input$col_POS,]
    }
    
    if (input$col_REF!= "All") {
      filtered_data <- filtered_data[filtered_data$REF == input$col_REF,]
    }
    
    if (input$col_ALT!= "All") {
      filtered_data <- filtered_data[filtered_data$ALT == input$col_ALT,]
    }
    
    if (input$col_SYMBOL!= "All") {
      filtered_data <- filtered_data[filtered_data$SYMBOL == input$col_SYMBOL,]
    }
    
    if (input$col_HGVSc!= "All") {
      filtered_data <- filtered_data[filtered_data$HGVSc == input$col_HGVSc,]
    }
    
    if (input$col_gnomAD_AF!= "All") {
      filtered_data <- filtered_data[filtered_data$gnomAD_AF == input$col_gnomAD_AF,]
    }
    
    if (input$col_SpliceAI_DS_AG!= "All") {
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_AG == input$col_SpliceAI_DS_AG,]
    }
    
    if (input$col_SpliceAI_DS_AL!= "All") {
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_AL == input$col_SpliceAI_DS_AL,]
    }
    
    if (input$col_SpliceAI_DS_DG!= "All") {
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_DG == input$col_SpliceAI_DS_DG,]
    }
    
    if (input$col_SpliceAI_DS_DL!= "All") {
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_DL == input$col_SpliceAI_DS_DL,]
    }
    
    if (input$col_mmsplice_delta_logit_psi!= "All") {
      filtered_data <- filtered_data[filtered_data$mmsplice_delta_logit_psi == input$col_mmsplice_delta_logit_psi,]
    }
    
    if (input$col_MaxEntScan_alt!= "All") {
      filtered_data <- filtered_data[filtered_data$MaxEntScan_alt == input$col_MaxEntScan_alt,]
    }
    
    if (input$col_MaxEntScan_diff!= "All") {
      filtered_data <- filtered_data[filtered_data$MaxEntScan_diff == input$col_MaxEntScan_diff,]
    }
    
    if (input$col_SQUIRLS!= "All") {
      filtered_data <- filtered_data[filtered_data$SQUIRLS == input$col_SQUIRLS,]
    }
    
    if (input$col_Type!= "All") {
      filtered_data <- filtered_data[filtered_data$Type == input$col_Type,]
    }
    
    if (input$col_QUAL!= "All") {
      filtered_data <- filtered_data[filtered_data$QUAL == input$col_QUAL,]
    }
    
    ### Show table of filtered data ###
    output$splice_table <- DT::renderDataTable(DT::datatable({
      filtered_data
    }))
    
  })
}

##### Run app #####
shinyApp(ui = ui, server = server)

##### Notes for Deployment #####
# library(rsconnect)
# rsconnect::deployApp('/data/jess_tmp/fh/Rshiny/splice') 
