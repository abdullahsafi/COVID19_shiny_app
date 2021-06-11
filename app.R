library(shiny)
library(shinydashboard)
library(ggplot2)
library(limma)
library(DT)
library(tidyverse)
library(plotly)
library(vroom)
library(dplyr)
library(stringr)
library(tools)

library(HiDimDA)
library(e1071)
library(randomForest)
load("train.Rdata")
load("train2.RData") #load training datasets and correspoding labels.

pdata_prevalidation = readRDS("pdata_prevalidation.RDS")
mdata_prevalidation = readRDS("mdata_prevalidation.RDS")
rdata_prevalidation = readRDS("rdata_prevalidation.RDS")

load("final_model.RData")
load("logistic_models.RData")

scaleMinMax <- function(x){
  (x - min(x)) / (max(x) - min(x))
}
#set up modelling
py = Combined_train_pdata_label
my = Combined_train_mdata_label 
ry = Combined_train_rdata_label 

pX <- t(Combined_train_pdata)
mX <- t(Combined_train_mdata)
rX <- t(Combined_train_rdata) # features in columns now



design_p <- model.matrix(~ factor(py))
fit_p <- lmFit(t(pX), design_p)
fit2_p <- eBayes(fit_p)
static_tstats_tT_p = topTable(fit2_p,coef = 2,number = Inf,sort.by = "t")

design_m <- model.matrix(~ factor(my))
fit_m <- lmFit(t(mX), design_m)
fit2_m <- eBayes(fit_m)
static_tstats_tT_m = topTable(fit2_m,coef = 2,number = Inf,sort.by = "t")

design_r <- model.matrix(~ factor(ry))
fit_r <- lmFit(t(rX), design_r)
fit2_r <- eBayes(fit_r)
static_tstats_tT_r = topTable(fit2_r,coef = 2,number = Inf,sort.by = "t")

pdata_samples = rownames(pX)
mdata_samples = rownames(mX)
rdata_samples = rownames(rX)

proteins = colnames(pX)
metabolomics = colnames(mX)
transcriptomics = colnames(rX)

#example of preparing fs on common_sample to rds
#source("model_topn_cv_visualise.R")
#source("model_topn_cv_visualise_recall.R")
#g_df = model_topn_cv_visualise_recall(topTable(fit2_g, coef=2,number = Inf, sort.by ="t"),gX,stack_label)
#saveRDS(g_df,"g_recall.rds")

## load training model



#top_features <- rownames(topTable(fit2_t, coef=2, sort.by ="t"))[1:50]
final_prob = c(0.5,0.5)
num = c(5,10,20,50,100)



ui <- dashboardPage(
  
  # Application title
  dashboardHeader(title = "Covid-19 risk calculator"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName = "0", icon = icon("table")),
      menuItem("Design", tabName = "1", icon = icon("chart-line")),
      menuItem("Prediction", tabName = "2", icon = icon("user"))
      
    )
  ),
  dashboardBody(
    tags$head(tags$style(
      HTML('.wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}')
    )),
    tabItems(
      tabItem(
        tabName = "0",
        fluidRow(
          tabBox(
            title = NULL, width  =12,
            id = "tabset0",height = "1500px",
            
            tabPanel("Overview", 
                     box(width = NULL, title = "Data Source", status = "primary",
                         uiOutput("ui1")
                     ),
                     
                     box(width = NULL, title = "Introduction", status = "primary", fluidRow(column(width = 12,
                                         div(style = "font-size:15px;",
                                             
                                         p('We have developed an omics-based risk calculator that enables the early identification of Covid-19 progression and therefore allows for early intervention. 
                                                                         The model classifies patients into “moderate” or “severe” Covid-19 risks. The features of this application include:'),
                                         tags$div(tags$ul(
                                           tags$li("Containing a multi-omic risk calculator."),
                                           tags$li("Deep insights into proteomic, metabolic, and transcriptomics platforms and their relationship Covid-19."),
                                           tags$li("Risk calculator uses a separate final model for each possible combination of omics platforms."),  
                                           tags$li("Our model and calculator can be adapted to users uploaded omics data")), style = "font-size: 15px")
                                         ))
                     )),box(width = NULL, title = "Disclaimer", status = "primary",
                            p("The interpretation of the results of the calculator is not recommend for those without appropriate medical and clinical knowledge. We do not accept any liability, including for any loss or damage, resulting from the reliance on the content, or for its accuracy and completeness.")
                            )
                     ),
            
            tabPanel("Proteomics",
                     fluidRow(column(width = 4,
                                     box(title = "Select Inputs", status = "primary",width = NULL,
                                         sliderInput(inputId = "variable_2", "Number of proteins selected: ", 1, length(proteins), 50)
                                     )
                     ),
                     column(width = 8,
                            box(width = NULL, title = "Top table ranked by moderate t-statistics", status = "primary",
                                DT::dataTableOutput("toptable_p")
                            )
                     )),
                     fluidRow(
                       column(width = 4,
                              box(title = "Select Protein:", status = "primary",width = NULL,
                                  selectInput(
                                    inputId = "protein",
                                    "Select protein:",
                                    choices = proteins,
                                    multiple = FALSE
                                  ))
                       ),
                       column(width = 8,
                              box(width = NULL, title = "Differential Expression", status = "primary",
                                  plotlyOutput("proteinDEPlot")
                              )
                       )
                     ),
                     fluidRow(
                       column(width = 4,
                              box(title = "Select a set of proteins:", status = "primary",width = NULL,
                                  selectInput(
                                    inputId = "protein_pca",
                                    "Select proteins:",
                                    choices = proteins,
                                    selected = c("P80098","P05231","P08727"),
                                    multiple = TRUE
                                  ))
                       ),
                       column(width = 8,
                              box(width = NULL, title = "PCA plot on the selected proteins", status = "primary",
                                  plotlyOutput("PCA_pPlot")
                              )
                       )
                     )
                     
            ),
            tabPanel("Metabolomics",
                     fluidRow(column(width = 4,
                                     box(title = "Select Inputs", status = "primary",width = NULL,
                                         sliderInput(inputId = "variable_3", "Number of features selected: ", 1, length(metabolomics), 50)
                                     )
                     ),
                     column(width = 8,
                            box(width = NULL, title = "Top table ranked by moderate t-statistics", status = "primary",
                                DT::dataTableOutput("toptable_m")
                            )
                     )),
                     fluidRow(
                       column(width = 4,
                              box(title = "Select metabolomic:", status = "primary",width = NULL,
                                  selectInput(
                                    inputId = "metabolomic",
                                    "Select metabolomic:",
                                    choices = metabolomics,
                                    selected = "1-(1-enyl-palmitoyl)-2-arachidonoyl-GPC (P-16:0/20:4)*",
                                    multiple = FALSE
                                  ))
                       ),
                       column(width = 8,
                              box(width = NULL, title = "Differential Expression", status = "primary",
                                  plotlyOutput("metabolomicDEPlot")
                              )
                       )
                     ),
                     fluidRow(
                       column(width = 4,
                              box(title = "Select a set of metabolites:", status = "primary",width = NULL,
                                  selectInput(
                                    inputId = "metabolite_pca",
                                    "Select metabolites:",
                                    choices = metabolomics,
                                    selected = c("S-methylcysteine sulfoxide","S-methylcysteine"),
                                    multiple = TRUE
                                  ))
                       ),
                       column(width = 8,
                              box(width = NULL, title = "PCA plot on selected metabolites", status = "primary",
                                  plotlyOutput("PCA_mPlot")
                              )
                       )
                     )
            ),
            tabPanel("Transcriptomics",
                     fluidRow(column(width = 4,
                                     box(title = "Select Inputs", status = "primary",width = NULL,
                                         sliderInput(inputId = "variable_0", "Number of features selected: ", 1, 200, 50)
                                     )
                     ),
                     column(width = 8,
                            box(width = NULL, title = "Top table ranked by moderate t-statistics", status = "primary",
                                DT::dataTableOutput("toptable_r")
                            )
                     )),
                     fluidRow(
                       column(width = 4,
                              box(title = "Select gene:", status = "primary",width = NULL,
                                  selectInput(
                                    inputId = "genomic",
                                    "Select gene:",
                                    choices = transcriptomics,
                                    #selected = "1-(1-enyl-palmitoyl)-2-arachidonoyl-GPC (P-16:0/20:4)*",
                                    multiple = FALSE
                                  ))
                       ),
                       column(width = 8,
                              box(width = NULL, title = "Differential Expression", status = "primary",
                                  plotlyOutput("genomicDEPlot")
                              )
                       )
                     ),
                     fluidRow(
                       column(width = 4,
                              box(title = "Select a set of transcriptomes:", status = "primary",width = NULL,
                                  selectInput(
                                    inputId = "transcriptome_pca",
                                    "Select transcriptomes:",
                                    choices = transcriptomics,
                                    selected = c("GRB10","UPP1"),
                                    multiple = TRUE
                                  ))
                       ),
                       column(width = 8,
                              box(width = NULL, title = "PCA plot on selected genes", status = "primary",
                                  plotlyOutput("PCA_rPlot")
                              )
                       )
                     )
            )
          )
        )
      ),
      tabItem(
        tabName = "1",
        fluidRow(
          tabBox(
            title = NULL, width = 12,
            # The id lets us use input$tabset1 on the server to find the current tab
            id = "tabset1", height = "1500px",
            
            # tabPanel("Single-platform Models",
            #          fluidRow(column(width = 12,
            #                          div(style = "font-size:15px;",
            #                           p('On each omics platform, we train models to distinguish between moderate/severe outcomes. 
            #                           To select the best, 5-fold cv (over 50 repetitions) was used.Further, to balance between bias/variance for each model, the number of features selected (ranked by t-statistic) was also varied.'),
            #                           p("We select models with good accuracy, but more importantly," ,em("high"),em(a("recall", href = "https://en.wikipedia.org/wiki/Precision_and_recall",target="_blank"))," - we are more concerned with detecting severe cases, as these are the more life-threatening.
            #                             "),
            #                           p("For each platform, the resultant cv recall and accuracy are shown below. Hence we choose:"),
            #                           p(code("Proteomics:"),"non-linear SVM, n = 40"),
            #                           p(code("Metalobomics:"),"non-linear SVM, n = 40"),
            #                           p(code("Genomics:"),"non-linear SVM, n =50")
            #         ))
            #         ),
            #         fluidRow(
            #               tabBox(id = "p_fs",title = "Performance on Proteomics",width = 12,height = "375px",
            #                      tabPanel("Recall",plotlyOutput("p_recall")),
            #                      tabPanel("Accuracy",plotlyOutput("p_acc")))
            #          
            #         ),
            #         fluidRow(
            #           tabBox(id = "m_fs",title = "Performance on Metalobomics",width = 12,height = "375px",
            #                  tabPanel("Recall",plotlyOutput("m_recall")),
            #                  tabPanel("Accuracy",plotlyOutput("m_acc")))
            #           
            #         ),
            #         fluidRow(
            #           tabBox(id = "g_fs",title = "Performance on Genomics",width = 12,height = "375px",
            #                  tabPanel("Recall",plotlyOutput("g_recall")),
            #                  tabPanel("Accuracy",plotlyOutput("g_acc")))
            #           
            #         )
            #        
            # ),
            tabPanel("Model Design", 
                     box(width = NULL,title = "Model Design Diagram",state = "primary",img(src='fig1screenshot.png', align = "left")))
          #   ,tabPanel("Final Model",
          #            fluidRow(column(width = 12,
          #                            div(style = "font-size:15px;",
          #                                p('On each omics platform, we choose the best model based on both recall and accuracy and calculate a corresponding prevalidation vector. The three prevalidation vectors is combined with some clinical variables to train a final model using KNN.')
          #                            )
          #            )
          #            ),
          #            fluidRow(column(width = 6,
          #                            box(title = "Select Input", status = "primary",width = NULL,
          #                                selectInput(
          #                                  inputId = "sample",
          #                                  "Select a sample:",
          #                                  choices = common_sample,
          #                                  multiple = FALSE
          #                                )
          #                            ),
          #                            valueBoxOutput("actual_label",width = 12)
          #                            
          #            ),
          #            column(width = 6,
          #             valueBoxOutput("final_prediction",width = 12),
          #             box(flexdashboard::gaugeOutput("percentage_chart"),width=12,
          #                       title="Risk Score"),
          #             
          #             valueBoxOutput("p_prediction",width = 12),
          #             valueBoxOutput("m_prediction",width = 12),
          #             valueBoxOutput("g_prediction",width = 12)
          #            )
          #            
          #         )
          #   
          # )
        )
      )
      ),
      tabItem(tabName = "2",
              tabBox(
                title = NULL, width  =12,
                id = "tabset0",height = "1500px", 
                tabPanel(title = "Prediction",fluid = TRUE,
                         fluidPage(
                           fluidRow(
                             box(width = NULL, title = "Instruction of uploading dataset", status = "primary", fluidRow(column(width = 12,
                                                                                                                                   div(style = "font-size:15px;",
                                                                                                                                       
                                                                                                                                       p('The current version of the shiny app only supports the datasets of following formats:'),
                                                                                                                                       tags$div(tags$ul(
                                                                                                                                         tags$li("The clinical data contains 3 columns with column names: sample_id,sex,age. This dataset is compulsory."),
                                                                                                                                         tags$li("The omics data should contains the top 100 features in the corresponding top tables ranked by moderate t-statistics for each platform"),
                                                                                                                                         tags$li("You need to reload the Shiny app to make a second prediction.")), style = "font-size: 15px")
                                                                                                                                   ))
                             ))
                           ),
                           fileInput("filec", NULL, buttonLabel = "Upload clinical data",accept = c(".csv", ".tsv")),
                           tableOutput("filecdetails"),
                           
                           fluidRow(
                             column(width = 4, 
                                    fileInput("filep", NULL, buttonLabel = "Proteomics: upload csv/tsv",accept = c(".csv", ".tsv",".rds")),
                                    #checkboxInput("filepcols", "Are features in columns?",value = FALSE),
                                    tableOutput("filepdetails")),
                             column(width = 4,
                                    fileInput("filem", NULL, buttonLabel = "Metalobomics: upload csv/tsv",accept = c(".csv", ".tsv")),
                                    #checkboxInput("filemcols", "Are features in columns?",value = FALSE),
                                    tableOutput("filemdetails")),
                             column(width = 4,
                                    fileInput("fileg", NULL, buttonLabel = "Transcriptomics: upload csv/tsv",accept = c(".csv", ".tsv")),
                                    #checkboxInput("filegcols", "Are features in columns?",value = FALSE),
                                    tableOutput("filegdetails"))
                           ),
                           
                           actionButton("predict", "Make prediction"),
                           hr(),
                           box(width = 12, title = "Predicted outcome", status = "primary",
                           DT::dataTableOutput("contents"))
                         )             
                )
                )
      )
    )
    
  )
)

server <- function(input, output) {
  
  tT_p = reactive({
    tT_p = round(topTable(fit2_p, coef=2, number = input$variable_2, sort.by ="t"),2)
    names = rownames(tT_p)
    #new_names = gsub(" ", "+", names)
    new_names = paste(paste(paste("<a href=\"https://www.uniprot.org/uniprot/", names,sep = ""),"\">",sep = ''),names,sep = "")
    new_names = paste(new_names,"</a>",sep = "")
    #'<a href="http://rstudio.com">RStudio</a>',
    rownames(tT_p) = new_names
    return(tT_p)
  })
  
  tT_m = reactive({
    tT_m = round(topTable(fit2_m, coef=2, number = input$variable_3, sort.by ="t"),2)
    names = rownames(tT_m)
    new_names = paste(paste(paste("<a href=\"https://hmdb.ca/unearth/q?utf8=%E2%9C%93&query=", names,sep = ""),"&searcher=metabolites&button=\">",sep = ''),names,sep = "")
    new_names = paste(new_names,"</a>",sep = "")
    #https://hmdb.ca/unearth/q?utf8=%E2%9C%93&query=1-stearoyl-2-oleoyl-GPE%20(18:0/18:1)&searcher=metabolites&button=
    rownames(tT_m) = new_names
    return(tT_m)
  })
  
  tT_r = reactive({
    tT_r = round(topTable(fit2_r, coef=2, number = input$variable_0, sort.by ="t"),2)
    return(tT_r)
  })
  
  output$PCA_pPlot <- renderPlotly({
    pca_protein = prcomp(pX[,input$protein_pca])
    df_pca <- data.frame(pca_protein$x)
    df_pca$y <- py
    # df_pca$y <- factor(y)
    p = ggplot(df_pca, aes(x = PC1, y = PC2, color = y)) + geom_point() + 
      theme_bw()
    ggplotly(p)
  })
  
  output$PCA_mPlot <- renderPlotly({
    pca_metabolite = prcomp(mX[,input$metabolite_pca])
    df_pca <- data.frame(pca_metabolite$x)
    df_pca$y <- my
    # df_pca$y <- factor(y)
    p = ggplot(df_pca, aes(x = PC1, y = PC2, color = y)) + geom_point() + 
      theme_bw()
    ggplotly(p)
  })
  
  output$PCA_rPlot <- renderPlotly({
    pca_transcriptome = prcomp(rX[,input$transcriptome_pca])
    df_pca <- data.frame(pca_transcriptome$x)
    df_pca$y <- ry
    # df_pca$y <- factor(y)
    p = ggplot(df_pca, aes(x = PC1, y = PC2, color = y)) + geom_point() + 
      theme_bw()
    ggplotly(p)
  })
  
  
  
  output$ui1 <- renderUI({
    p("To train our models, we use data from the following paper by Su et. al: ",em("Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19"), " ( 2020, Cell 183, 1479–1495. DOI: ",
a("https://doi.org/10.1016/j.cell.2020.10.037", href = "https://doi.org/10.1016/j.cell.2020.10.037",target="_blank"),"), which contains 3 platforms, proteomics, metabolomics, and transcriptomics. We supplement each of these platforms using the other three papers. ")
  })
  
  output$toptable_p <- DT::renderDataTable(datatable(tT_p(), options = list(pageLength = 5),escape = FALSE))
  
  output$proteinDEPlot <- renderPlotly({
    expression <- data.frame(protein = pX[,input$protein],
                             patient = py)
    colnames(expression) <- c("protein", "groups")
    mtitle <- paste("Expression of", input$protein)
    
    p <- ggplot(expression, aes(x = groups, y = protein, colour = groups)) + 
      geom_boxplot() +
      labs(title = mtitle, x = "patient outcome", y = "log2 expression") + 
      theme_minimal() + theme(legend.position = "none")
                                                                                 
      ggplotly(p) %>% layout(title = list(size = 1))
  })
  
  output$toptable_m <- DT::renderDataTable(datatable(tT_m(), options = list(pageLength = 5),escape = FALSE))
  
  output$metabolomicDEPlot <- renderPlotly({
    expression <- data.frame(metabolomic = mX[,input$metabolomic],
                             patient = my)
    colnames(expression) <- c("metabolomic", "groups")
    mtitle <- stringr::str_wrap(paste("Expression of", input$metabolomic),width = 60)
    
    m <- ggplot(expression, aes(x = groups, y = metabolomic, colour = groups)) + 
      geom_boxplot() +
      labs(title = mtitle, x = "patient outcome", y = "log2 expression")   + 
      theme_minimal() + theme(legend.position = "none") 
    
    ggplotly(m)
  })
  
  output$toptable_r <- DT::renderDataTable(datatable(tT_r()[1:1000,], options = list(pageLength = 5)))
  
  output$genomicDEPlot <- renderPlotly({
    expression <- data.frame(genomic = rX[,input$genomic],
                             patient = ry)
    colnames(expression) <- c("genomic", "groups")
    mtitle <- stringr::str_wrap(paste("Expression of", input$genomic),width = 60)
    
    g <- ggplot(expression, aes(x = groups, y = genomic, colour = groups)) + 
      geom_boxplot() +
      labs(title = mtitle, x = "patient outcome", y = "log2 expression")   + 
      theme_minimal() + theme(legend.position = "none") 
    
    ggplotly(g)
  })
  
  # output$actual_label <- renderValueBox({ ###update
  #   # The following code runs inside the database.
  #   # pull() bring the results into R, which then
  #   # it's piped directly to a valueBox()
  #   result = stack_label[which(input$sample == common_sample)]
  #   valueBox(result, icon = icon("list"),
  #            color = "green",
  #            subtitle = "Actual Label",
  #            width = 12
  #   )
  # })
  
  # output$final_prediction <- renderValueBox({ ###update
  #   # The following code runs inside the database.
  #   # pull() bring the results into R, which then
  #   # it's piped directly to a valueBox()
  #   result = ifelse(predicted[which(input$sample == common_sample)] == 0, "moderate","severe")
  #   valueBox(result, icon = icon("list"),
  #            color = "purple",
  #            subtitle = "Final Prediction",
  #            width = 12
  #   )
  # })
  # 
  # output$m_prediction <- renderValueBox({ ###update
  #   # The following code runs inside the database.
  #   # pull() bring the results into R, which then
  #   # it's piped directly to a valueBox()
  #   result = prevalidation_df[which(input$sample == common_sample),4]
  #   valueBox(paste0(round(result*100,0),"%"), icon = icon("list"),
  #            color = "blue",
  #            subtitle = "to be severe case predicted from Metabolomics",
  #            width = 12
  #   )
  # })
  # 
  # output$p_prediction <- renderValueBox({ ###update
  #   # The following code runs inside the database.
  #   # pull() bring the results into R, which then
  #   # it's piped directly to a valueBox()
  #   result = prevalidation_df[which(input$sample == common_sample),1]
  #   valueBox(paste0(round(result*100,0),"%"), icon = icon("list"),
  #            color = "yellow",
  #            subtitle = "to be severe case predicted from Proteomics",
  #            width = 12
  #   )
  # })
  # 
  # output$g_prediction <- renderValueBox({ ###update
  #   # The following code runs inside the database.
  #   # pull() bring the results into R, which then
  #   # it's piped directly to a valueBox()
  #   result = prevalidation_df[which(input$sample == common_sample),5]
  #   valueBox(paste0(round(result*100,0),"%"), icon = icon("list"),
  #            color = "red",
  #            subtitle = "to be severe case predicted from Genomics",
  #            width = 12
  #   )
  # })
  # 
  # output$percentage_chart <- flexdashboard::renderGauge({
  #   result = final_prob[which(input$sample == common_sample)]
  #   flexdashboard::gauge(round(as.numeric(result),2)*100, min = 0, max = 100, symbol = '%')
  # })
  
  
  ### Uploading
  data_c <- reactive({
    req(input$filec)
    
    ext <- tools::file_ext(input$filec$name)
    
    switch(ext,
           csv = vroom::vroom(input$filec$datapath, delim = ","),
           tsv = vroom::vroom(input$filec$datapath, delim = "\t"),
           rds = readRDS(input$filec$datapath),
           validate("Invalid file; Please upload a .csv, .tsv or .rds")
    )
  })
  
  data_p <- reactive({
    req(input$filep)
    
    ext <- tools::file_ext(input$filep$name)
    switch(ext,
           csv = vroom::vroom(input$filep$datapath, delim = ","),
           tsv = vroom::vroom(input$filep$datapath, delim = "\t"),
           rds = readRDS(input$filep$datapath),
           validate("Invalid file; Please upload a .csv, .tsv or .rds")
    )
  })
  
  data_m <- reactive({
    req(input$filem)
    
    ext <- tools::file_ext(input$filem$name)
    switch(ext,
           csv = vroom::vroom(input$filem$datapath, delim = ","),
           tsv = vroom::vroom(input$filem$datapath, delim = "\t"),
           rds = readRDS(input$filem$datapath),
           validate("Invalid file; Please upload a .csv, .tsv or .rds")
    )
  })
  
  data_g <- reactive({
    req(input$fileg)
    
    ext <- tools::file_ext(input$fileg$name)
    switch(ext,
           csv = vroom::vroom(input$fileg$datapath, delim = ","),
           tsv = vroom::vroom(input$fileg$datapath, delim = "\t"),
           rds = readRDS(input$fileg$datapath),
           validate("Invalid file; Please upload a .csv, .tsv or .rds")
    )
  })
  
  output$filecdetails <- renderTable(input$filec[1:3])
  output$filepdetails <- renderTable(input$filep[1:3])
  output$filemdetails <- renderTable(input$filem[1:3])
  output$filegdetails <- renderTable(input$fileg[1:3])
  
  
  testing_data <- eventReactive(input$predict, {
    
    if (is.null(input$filec)){
      validate("Please provide clinical data")
    }
    
    if (is.null(input$filep) & is.null(input$filem) & is.null(input$fileg)){
      validate("Please provide one of proteomics, metalobomics or genomics")
    }
    
    # at this point, clinical is available, and one of the platforms are available
    
    #load clinical
    clinical_test = data_c()
    clinical_test = clinical_test %>% dplyr::select(sample_id,sex,age)
    
    platforms_predictions = NULL
    #if protein 
    if  (!is.null(input$filep)){
      new_data = as.data.frame(data_p())
      
      #for now, only rds
      # if (input$filepcols  == FALSE){ #if features aren't in cols
      #   new_data = t(new_data)
      # }
      #adjustments
      
      #order
      new_data <- new_data[clinical_test$sample_id, ]
      
      #new_data = as.data.frame(apply(new_data, 2, scaleMinMax))
      #for proteomics: non-linear SVM (n = 40), SVM (n = 50), RF (n = 20)
      svm_data = new_data %>% dplyr::select(rownames(static_tstats_tT_p)[1:40])
      nl_svm_data = new_data %>% dplyr::select(rownames(static_tstats_tT_p)[1:50])
      rf_data = new_data %>% dplyr::select(rownames(static_tstats_tT_p)[1:20])
      
      pdata_predicted_nl_svm <- predict(pdata_trained_nl_svm, nl_svm_data)
      pdata_predicted_svm <- predict(pdata_trained_svm, svm_data)
      pdata_predicted_rf <- predict(pdata_trained_rf, rf_data)
      pdata_predicted_df = data.frame(cbind(pdata_predicted_nl_svm,pdata_predicted_svm,pdata_predicted_rf))
      #majority voting
      pdata_predicted_df <- sapply(pdata_predicted_df, as.numeric)
      pdata_voted = ifelse(rowMeans(pdata_predicted_df) > 1.5, "severe","moderate")
      # predicted_prob = attr(predicted_nl_svm,"probabilities")[,"severe"]
      # 
       if (is.null(platforms_predictions)){
         platforms_predictions = data.frame(p_preval = pdata_voted)
       }else{
         platforms_predictions = cbind(platforms_predictions,p_preval = pdata_voted)
       }
      
    }
    
    #if metalobomics
    if  (! is.null(input$filem)){
      new_data = as.data.frame(data_m())
      
      #for now, only rds
      #if (input$filemcols  == FALSE){ #if features aren't in cols
      #  new_data = t(new_data)
      #}
      #adjustments
      #order
      new_data <- new_data[clinical_test$sample_id,]
      #for metabolomics: non-linear SVM (n = 50), SVM (n = 50), DLDA (n = 10)
      svm_data = new_data %>% dplyr::select(rownames(static_tstats_tT_m)[1:50])
      nl_svm_data = new_data %>% dplyr::select(rownames(static_tstats_tT_m)[1:50])
      dlda_data = new_data %>% dplyr::select(rownames(static_tstats_tT_m)[1:10])
      
      mdata_predicted_nl_svm <- predict(mdata_trained_nl_svm, nl_svm_data)
      mdata_predicted_svm <- predict(mdata_trained_svm, svm_data)
      mdata_predicted_dlda <- predict(mdata_trained_dlda, dlda_data,grpcodes=levels(factor(Combined_train_mdata_label)))$class
      mdata_predicted_df = data.frame(cbind(mdata_predicted_nl_svm,mdata_predicted_svm,mdata_predicted_dlda))
      #majority voting
      mdata_predicted_df <- sapply(mdata_predicted_df, as.numeric)
      mdata_voted = ifelse(rowMeans(mdata_predicted_df) > 1.5, "severe","moderate")
      
      if (is.null(platforms_predictions)){
        platforms_predictions = data.frame(m_preval = mdata_voted) #define col name
      } else {
        platforms_predictions = cbind(platforms_predictions,m_preval = mdata_voted)
      }
    }
    
    
    ###
    #if genomics
    if  (! is.null(input$fileg)){
      new_data = as.data.frame(data_g())
      
      #for now, only rds
      # if (input$filegcols  == FALSE){ #if features aren't in cols
      #   new_data = t(new_data)
      # }
      #adjustments
      #order
      new_data <- new_data[clinical_test$sample_id, ]
      
      #for transcriptomics, nl-svm, svm and rf (n = 50)
      svm_data = new_data %>% dplyr::select(rownames(static_tstats_tT_r)[1:50])
      nl_svm_data = new_data %>% dplyr::select(rownames(static_tstats_tT_r)[1:50])
      rf_data = new_data %>% dplyr::select(rownames(static_tstats_tT_r)[1:50])
      
      rdata_predicted_nl_svm <- predict(rdata_trained_nl_svm, nl_svm_data)
      rdata_predicted_svm <- predict(rdata_trained_svm, svm_data)
      rdata_predicted_rf <- predict(rdata_trained_rf, rf_data)
      rdata_predicted_df = data.frame(cbind(rdata_predicted_nl_svm,rdata_predicted_svm,rdata_predicted_rf))
      #majority voting
      rdata_predicted_df <- sapply(rdata_predicted_df, as.numeric)
      rdata_voted = ifelse(rowMeans(rdata_predicted_df) > 1.5, "severe","moderate")
      
      if (is.null(platforms_predictions)){
        platforms_predictions = data.frame(r_preval = rdata_voted)
      }else{
        platforms_predictions = cbind(platforms_predictions,r_preval = rdata_voted)
      }
    }
    #return(clinical_test)
    return(cbind(clinical_test,platforms_predictions))
  })
  
  
  
  output$contents <- DT::renderDataTable({
    if (!is.null(input$filep) & is.null(input$filem) & is.null(input$fileg)){#only proteomics
      predicted_glm <-  predict(glm_p, testing_data(), type = "response")
      predicted_label <- ifelse (predicted_glm > 0.5 ,"severe", "moderate")
      result = cbind(testing_data() %>% dplyr::select(sample_id,sex,age),
                     final_prediction = predicted_label,
                     model_probability = predicted_glm)
      return(DT::datatable(
        result %>% 
          mutate_if(is.numeric, round,digits=1),
        options = list(pageLength = 5,
                       autoWidth = FALSE, scrollX = TRUE
        )))
      #df = testing_data()
      #return(cbind(testing_data(),final_prediction = if_else(df$p_pred_nl_svm > 0.5, "severe", "moderate")))
    }
    
    if (is.null(input$filep) & !is.null(input$filem) & is.null(input$fileg)){#only metalobomics
      predicted_glm <-  predict(glm_m, testing_data(), type = "response")
      predicted_label <- ifelse (predicted_glm > 0.5 ,"severe", "moderate")
      
      result = cbind(testing_data() %>% dplyr::select(sample_id,sex,age),
                     final_prediction = predicted_label,
                     model_probability = predicted_glm)
      return(DT::datatable(
        result %>% 
          mutate_if(is.numeric, round,digits=1),
        options = list(pageLength = 5,
                       autoWidth = FALSE, scrollX = TRUE
        )))
      #df = testing_data()
      #return(cbind(testing_data(),final_prediction = if_else(df$m_pred_nl_svm > 0.5, "severe", "moderate")))
    }
    
    if (!is.null(input$filep) & !is.null(input$filem) & is.null(input$fileg)){#proteomics +  metalobomics
      predicted_glm <-  predict(glm_pm, testing_data(), type = "response")
      predicted_label <- ifelse (predicted_glm > 0.5 ,"severe", "moderate")
      result = cbind(testing_data() %>% dplyr::select(sample_id,sex,age),
                     final_prediction = predicted_label,
                     model_probability = predicted_glm)
      return(DT::datatable(
        result %>% 
          mutate_if(is.numeric, round,digits=1),
        options = list(pageLength = 5,
                       autoWidth = FALSE, scrollX = TRUE
        )))
      #df = testing_data()
      #vote = rowMeans(df %>% dplyr::select(p_pred_nl_svm,m_pred_nl_svm))
      #return(cbind(testing_data(),final_prediction = if_else(vote > 0.5, "severe", "moderate")))
    }
    
    if (is.null(input$filep) & is.null(input$filem) & !is.null(input$fileg)){#genomics
      predicted_glm <-  predict(glm_r, testing_data(), type = "response")
      predicted_label <- ifelse (predicted_glm > 0.5 ,"severe", "moderate")
      result = cbind(testing_data() %>% dplyr::select(sample_id,sex,age),
                     final_prediction = predicted_label,
                     model_probability = predicted_glm)
      return(DT::datatable(
        result %>% 
          mutate_if(is.numeric, round,digits=1),
        options = list(pageLength = 5,
                       autoWidth = FALSE, scrollX = TRUE
        )))
      #df = testing_data()
      #return(cbind(testing_data(),final_prediction = if_else(df$g_pred_nl_svm > 0.5, "severe", "moderate")))
    }
    
    if (!is.null(input$filep) & is.null(input$filem) & !is.null(input$fileg)){#genomics + proteomics
      predicted_glm <-  predict(glm_pr, testing_data(), type = "response")
      predicted_label <- ifelse (predicted_glm > 0.5 ,"severe", "moderate")
      result = cbind(testing_data() %>% dplyr::select(sample_id,sex,age),
                     final_prediction = predicted_label,
                     model_probability = predicted_glm)
      return(DT::datatable(
        result %>% 
          mutate_if(is.numeric, round,digits=1),
        options = list(pageLength = 5,
                       autoWidth = FALSE, scrollX = TRUE
        )))
      #df = testing_data()
      #vote = rowMeans(df %>% dplyr::select(p_pred_nl_svm,g_pred_nl_svm))
      #return(cbind(testing_data(),final_prediction = if_else(vote > 0.5, "severe", "moderate")))
    }
    
    if (is.null(input$filep) & !is.null(input$filem) & !is.null(input$fileg)){#genomics + metalobomics
      predicted_glm <-  predict(glm_mr, testing_data(), type = "response")
      predicted_label <- ifelse (predicted_glm > 0.5 ,"severe", "moderate")
      result = cbind(testing_data() %>% dplyr::select(sample_id,sex,age),
                     final_prediction = predicted_label,
                     model_probability = predicted_glm)
      return(DT::datatable(
        result %>% 
          mutate_if(is.numeric, round,digits=1),
        options = list(pageLength = 5,
                       autoWidth = FALSE, scrollX = TRUE
        )))
      #df = testing_data()
      #vote = rowMeans(df %>% dplyr::select(g_pred_nl_svm,m_pred_nl_svm))
      #return(cbind(testing_data(),final_prediction = if_else(vote > 0.5, "severe", "moderate")))
    }
    
    if (!is.null(input$filep) & !is.null(input$filem) & !is.null(input$fileg)){#all in
      predicted_glm <-  predict(glm_pmr, testing_data(), type = "response")
      predicted_label <- ifelse (predicted_glm > 0.5 ,"severe", "moderate")
      result = cbind(testing_data() %>% dplyr::select(sample_id,sex,age),
                     final_prediction = predicted_label,
                     model_probability = predicted_glm)
      return(DT::datatable(
        result %>% 
          mutate_if(is.numeric, round,digits=1),
        options = list(pageLength = 5,
                       autoWidth = FALSE, scrollX = TRUE
                       )))
      #df = testing_data()
      #vote = rowMeans(df %>% dplyr::select(p_pred_nl_svm,m_pred_nl_svm,g_pred_nl_svm))
      #return(cbind(testing_data(),final_prediction = if_else(vote > 0.5, "severe", "moderate")))
    }
  })
  
  
  
}

shinyApp(ui, server)