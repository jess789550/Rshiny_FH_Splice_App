# DEPLOYED VERSION #

##### Load libraries ######
library(shiny)
library(DT)
library(shinydashboard) # for box()
library(ggplot2)
library(stringr)

##### Define UI: https://shiny.posit.co/r/gallery/widgets/basic-datatable/ #####
ui <- dashboardPage(
  ### Title ##
  dashboardHeader(
    title="Exome Splice Site Prediction Results"
  ),
  
  ### Create a sidebar panel in the UI for adjusting parameters ###
  dashboardSidebar(
    
    # Allow user to upload worklist
    # fileInput("file1", "Choose CSV File", accept = ".csv"),
    # checkboxInput("header", "Header", TRUE),
    
    # Allow user to select worklist
    selectInput("worklist", "Choose a worklist:",
                choices = c("1234567", "1234568", "1234569", "1234570", "1234571",
                            "1234572", "1234573", "1234574")),
    
    # Choose SpliceAI cutoff
    sliderInput("SpliceAI", "SpliceAI cutoff:", min = 0, max = 1, value = 0.2, step = 0.1),
    
    # Choose MES cutoff
    # sliderInput("MES", "MES cutoff:", min = 0, max = 10, value = 6.2, step = 0.1),
    selectInput("MES", "MES cutoff:", choices = c("None", "Low","High"), selected="Low"),
    
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
                                                                         "0.002","0.005", "0.01"), selected="0.01"),
    
    # Detected by GeneSplicer?
    radioButtons("GeneSplicer", label = "Must be detected by GeneSplicer?",
                 choices = c("Yes", "No"), selected = "No"),
    
    # Choose what to plot
    selectInput("plot", "FDR plot:", choices = c("SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", 
                                                 "SpliceAI_DS_DL", "MaxEntScan_alt", "SQUIRLS", "mmsplice_delta_logit_psi")),
    
    # Add button for user to press to initiate run
    actionButton(inputId = "Submit", label = "Submit"),
    
    # Add button to reset slider inputs
    actionButton(inputId = "Reset", label = "Reset")
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
                         fluidRow(box(width=12, style='width:1000px',
                                      fluidRow(
                                        column(2, selectizeInput(inputId="col_file_id", label="File:", choices="All")),
                                        column(2, selectizeInput(inputId="col_CHROM", label="Chromosome:", choices="All")),
                                        column(2,selectizeInput(inputId="col_POS", label="Position:", 
                                                                choices=c("All", "0-100000",  "100000-30100000",  "30100000-60100000",  "60100000-90100000", 
                                                                          "90100000-120100000", "120100000-150100000", "150100000-180100000", 
                                                                          "180100000-210100000", "210100000-240100000", "240100000-270100000"))),
                                        column(2,selectizeInput(inputId="col_REF", label="Reference allele:", choices="All")),
                                        column(2,selectizeInput(inputId="col_ALT", label="Alternative allele:", choices="All")),
                                        column(2,selectizeInput(inputId="col_SYMBOL", label="Gene symbol:", 
                                                                choices=c("All", "TTN", "MYBPC3", "LDLR", "MYH7", "BRCA2", 
                                                                          "APOB", "BRCA1", "FLNC", "SCN5A", "VWF"))), # Top 10 most freq
                                      ),
                                      fluidRow(
                                        # column(2, selectizeInput(inputId="col_HGVSc", label="HGVSc nomenclature:", choices="All")),
                                        # column(2, selectizeInput(inputId="col_gnomAD_AF", label="gnomAD allele frequency:", 
                                        #                          choices=c("All", 0-0.25, 0.25-0.5, 0.5-0.75, 0.75-1))),
                                        column(2, selectizeInput(inputId="col_SpliceAI_DS_AG", label="SpliceAI Acceptor Gain:", 
                                                                 choices=c("All", "0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                                                                           "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"))),
                                        column(2, selectizeInput(inputId="col_SpliceAI_DS_AL", label="SpliceAI Acceptor Loss:", 
                                                                 choices=c("All", "0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                                                                           "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"))),
                                        column(2, selectizeInput(inputId="col_SpliceAI_DS_DG", label="SpliceAI Donor Gain:", 
                                                                 choices=c("All", "0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                                                                           "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"))),
                                        column(2, selectizeInput(inputId="col_SpliceAI_DS_DL", label="SpliceAI Donor Loss:", 
                                                                 choices=c("All", "0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                                                                           "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"))),
                                        column(2, selectizeInput(inputId="col_mmsplice_delta_logit_psi", label="MMSplice:", 
                                                                 choices=c("All", "Less than 0", "0-0.5", "0.5-1", "1-1.5", "1.5-2"))),
                                      ),
                                      fluidRow(
                                        column(2, selectizeInput(inputId="col_MaxEntScan_alt", label="MES alt:", 
                                                                 choices=c("All", "Less than 0", "0-5", "5-10", "10-15", "15-20"))),
                                        column(2, selectizeInput(inputId="col_MaxEntScan_diff", label="MES diff:", 
                                                                 choices=c("All", "Less than 0", "0-5", "5-10", "10-15", "15-20"))),
                                        column(2, selectizeInput(inputId="col_SQUIRLS", label="SQUIRLS:", 
                                                                 choices=c("All", "0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                                                                           "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"))),
                                        column(2, selectizeInput(inputId="col_Type", label="True positive or false positive:", choices="All")),
                                        column(2, selectizeInput(inputId="col_QUAL", label="VCF Quality score:", 
                                                                 choices=c("All", "0-3000", "3000-6000", "6000-9000", "9000-12000", "12000-15000",
                                                                           "15000-18000", "18000-21000", "21000-24000", "24000-27000", "27000-30000"))),
                                      ),
                                      fluidRow(
                                        column(2, actionButton(inputId = "Filter", label = "Filter"))
                                      ))),
                         br(),
                         fluidRow(box(width=10, style='width:800px;overflow-x: scroll; overflow-y: scroll;',
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
    #file <- paste(input$worklist, ".csv", sep="")
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    #original_data <<- read.csv(file)  # need original_data downstream
    original_data <<- read.csv(file$datapath, header = input$header)
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
    
    # updateSelectizeInput(session, "col_gnomAD_AF",
    #                   choices = c("All", unique(as.character(data$gnomAD_AF))))
    
    # updateSelectizeInput(session, "col_SpliceAI_DS_AG",
    #                   choices = c("All", unique(as.character(data$SpliceAI_DS_AG))))
    # 
    # updateSelectizeInput(session, "col_SpliceAI_DS_AL",
    #                   choices = c("All", unique(as.character(data$SpliceAI_DS_AL))))
    # 
    # updateSelectizeInput(session, "col_SpliceAI_DS_DG",
    #                   choices = c("All", unique(as.character(data$SpliceAI_DS_DG))))
    # 
    # updateSelectizeInput(session, "col_SpliceAI_DS_DL",
    #                   choices = c("All", unique(as.character(data$SpliceAI_DS_DL))))
    
    # updateSelectizeInput(session, "col_mmsplice_delta_logit_psi",
    #                   choices = c("All", unique(as.character(data$mmsplice_delta_logit_psi))))
    
    # updateSelectizeInput(session, "col_MaxEntScan_alt",
    #                   choices = c("All", unique(as.character(data$MaxEntScan_alt))))
    # 
    # updateSelectizeInput(session, "col_MaxEntScan_diff",
    #                   choices = c("All", unique(as.character(data$MaxEntScan_diff))))
    
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
      null_freq_data <- data[which(data$gnomAD_AF=='-'),]
      dataFilteredgnomAD <- data[which(data$gnomAD_AF!='-'),]
      dataFilteredgnomAD$gnomAD_AF <- as.numeric(dataFilteredgnomAD$gnomAD_AF)
      data <- subset(dataFilteredgnomAD, gnomAD_AF <= freq)
      data <- rbind(data, null_freq_data)
    }
    
    # Filter GeneSplicer
    if (input$GeneSplicer == "Yes") {
      data <- subset(data, GeneSplicer_score != "-")
    }
    
    filtered_data <<- data  # set global variable so filtered_data can be accessed below
    
    ### Show table of filtered data ###
    output$splice_table <- DT::renderDataTable(DT::datatable({
      data
    }))
    
    ### Performance metrics TP, FP, FDR = FP / (FP + TP) ###
    output$metrics_table <- DT::renderDataTable(DT::datatable({
      
      # Get metrics
      TP <- nrow(filtered_data[filtered_data$Type == "TP",])
      FP <- nrow(filtered_data[filtered_data$Type == "FP",])
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
    plot_func <<- function(dataset, column) {
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
        plot_func(filtered_data, input$plot)
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
  
  ### Filter table of splice variants and prediction scores after clicking Filter button ###
  observeEvent(input$Filter, {
    # Filter by column values
    # Filter file names/sample IDs
    if (input$col_file_id != "All") {
      filtered_data <- filtered_data[filtered_data$file_id == input$col_file_id,]
    }
    
    # Filter chromosome
    if (input$col_CHROM!= "All") {
      filtered_data <- filtered_data[filtered_data$CHROM == input$col_CHROM,]
    }
    
    # Filter position
    if (input$col_POS!= "All") {
      sep <- str_split(input$col_POS, "-")
      start <- sep[[1]][1]
      end <- sep[[1]][2]
      filtered_data <- filtered_data[filtered_data$POS >= start,]
      filtered_data <- filtered_data[filtered_data$POS <= end,]
    }
    
    # Filter reference allele
    if (input$col_REF!= "All") {
      filtered_data <- filtered_data[filtered_data$REF == input$col_REF,]
    }
    
    # Filter alternative allele
    if (input$col_ALT!= "All") {
      filtered_data <- filtered_data[filtered_data$ALT == input$col_ALT,]
    }
    
    # Filter gene symbol
    if (input$col_SYMBOL!= "All") {
      filtered_data <- filtered_data[filtered_data$SYMBOL == input$col_SYMBOL,]
    }
    
    # Filter HGVSc transcript
    
    # if (input$col_HGVSc!= "All") {
    #   filtered_data <- filtered_data[filtered_data$HGVSc == input$col_HGVSc,]
    # }
    
    # Filter gnomad allele frequency column
    # if (input$col_gnomAD_AF!= "All") {
    #   filtered_data <- filtered_data[filtered_data$gnomAD_AF == input$col_gnomAD_AF,]
    # }
    
    # Filter SpliceAI columns
    if (input$col_SpliceAI_DS_AG!= "All") {
      sep <- str_split(input$col_SpliceAI_DS_AG, "-")
      start <- sep[[1]][1]
      end <- sep[[1]][2]
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_AG >= start,]
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_AG <= end,]
    }
    
    if (input$col_SpliceAI_DS_AL!= "All") {
      sep <- str_split(input$col_SpliceAI_DS_AL, "-")
      start <- sep[[1]][1]
      end <- sep[[1]][2]
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_AL >= start,]
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_AL <= end,]
    }
    
    if (input$col_SpliceAI_DS_DG!= "All") {
      sep <- str_split(input$col_SpliceAI_DS_DG, "-")
      start <- sep[[1]][1]
      end <- sep[[1]][2]
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_DG >= start,]
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_DG <= end,]
    }
    
    if (input$col_SpliceAI_DS_DL!= "All") {
      sep <- str_split(input$col_SpliceAI_DS_DL, "-")
      start <- sep[[1]][1]
      end <- sep[[1]][2]
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_DL >= start,]
      filtered_data <- filtered_data[filtered_data$SpliceAI_DS_DL <= end,]
    }
    
    # Filter MMSplice column
    if (input$col_mmsplice_delta_logit_psi!= "All") {
      if (input$col_mmsplice_delta_logit_psi == "Less than 0") {
        filtered_data <- filtered_data[filtered_data$mmsplice_delta_logit_psi <= 0,]
      } else {
        sep <- str_split(input$col_mmsplice_delta_logit_psi, "-")
        start <- sep[[1]][1]
        end <- sep[[1]][2]
        filtered_data <- filtered_data[filtered_data$mmsplice_delta_logit_psi >= start,]
        filtered_data <- filtered_data[filtered_data$mmsplice_delta_logit_psi <= end,]
      }
    }
    
    # Filter MES columns
    if (input$col_MaxEntScan_alt != "All") {
      if (input$col_MaxEntScan_alt == "Less than 0") {
        filtered_data <- filtered_data[filtered_data$MaxEntScan_alt <= 0,]
      } else {
        sep <- str_split(input$col_MaxEntScan_alt, "-")
        start <- sep[[1]][1]
        end <- sep[[1]][2]
        filtered_data <- filtered_data[filtered_data$MaxEntScan_alt >= start,]
        filtered_data <- filtered_data[filtered_data$MaxEntScan_alt <= end,]
      }
    }
    
    if (input$col_MaxEntScan_diff != "All") {
      if (input$col_MaxEntScan_diff == "Less than 0") {
        filtered_data <- filtered_data[filtered_data$MaxEntScan_diff <= 0,]
      } else {
        sep <- str_split(input$col_MaxEntScan_diff, "-")
        start <- sep[[1]][1]
        end <- sep[[1]][2]
        filtered_data <- filtered_data[filtered_data$MaxEntScan_diff >= start,]
        filtered_data <- filtered_data[filtered_data$MaxEntScan_diff <= end,]
      }
    }
    
    # Filter SQUIRLS column
    if (input$col_SQUIRLS!= "All") {
      sep <- str_split(input$col_SQUIRLS, "-")
      start <- sep[[1]][1]
      end <- sep[[1]][2]
      filtered_data <- filtered_data[filtered_data$SQUIRLS >= start,]
      filtered_data <- filtered_data[filtered_data$SQUIRLS <= end,]
    }
    
    # Filter TP/FP column
    if (input$col_Type!= "All") {
      filtered_data <- filtered_data[filtered_data$Type == input$col_Type,]
    }
    
    # Filter VCF quality score column
    if (input$col_QUAL!= "All") {
      sep <- str_split(input$col_QUAL, "-")
      start <- sep[[1]][1]
      end <- sep[[1]][2]
      filtered_data <- filtered_data[filtered_data$QUAL >= start,]
      filtered_data <- filtered_data[filtered_data$QUAL <= end,]
    }
    
    # Show table of filtered data #
    output$splice_table <- DT::renderDataTable(DT::datatable({
      filtered_data
    }))
    
    ### Refresh Performance metrics TP, FP, FDR = FP / (FP + TP) ###
    output$metrics_table <- DT::renderDataTable(DT::datatable({
      
      # Get metrics
      TP <- nrow(filtered_data[filtered_data$Type == "TP",])
      FP <- nrow(filtered_data[filtered_data$Type == "FP",])
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
    
    ### Refresh TP-FP trade-off plot ###
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
        plot_func(filtered_data, input$plot)
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
  
  ### Reset slider inputs ###
  observeEvent(input$Reset, {
    
    # Reset SpliceAI cutoff
    updateSliderInput(session, "SpliceAI", value = 0.3)
    
    # Reset MES cutoff
    updateSelectInput(session, "MES", selected = "Low")
    
    # Reset SQUIRLS cutoff
    updateSliderInput(session, "SQUIRLS", value = 0.5)
    
    # Reset MMSplice cutoff
    updateSliderInput(session, "MMSplice", value = 0.5)
    
    # Reset gnomAD cutoff
    updateSelectInput(session, "gnomAD", selected = "0.01")
    
    # Reset Detected by GeneSplicer?
    updateRadioButtons(session, "GeneSplicer", selected = "No")
  })
}

##### Run app #####
shinyApp(ui = ui, server = server)

##### Notes for Deployment #####
# library(rsconnect)
# rsconnect::deployApp('/data/jess_tmp/exome/Rshiny/exome_splice') 
