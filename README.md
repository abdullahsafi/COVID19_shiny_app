# Covid_C1_shinyapp
This app explores multi-omics dataset for severe/moderate covid-19 cases and allows the users to upload their own dataset for prediction.
Run the app locally in Rstudio on laptop by

``` r
library(shiny)
library(shinydashboard)
library(ggplot2)
library(limma)
library(DT)
library(tidyverse)
library(plotly)
library(vroom)
library(HiDimDA)
library(e1071)
library(randomForest)
library(dplyr)
library(stringr)
library(tools)

shiny::runGitHub(
  repo = "Covid_C1_shinyapp", 
  username = "YYAO3610", 
  ref = "main")
```

The demo_uploading_file folder contains 2 different demo uploading file using the test set from Su et al/Shen et al. You can use these files as testing uploading files. pdata,mdata,rdata refer to proteomics,metabolomics,transcriptomics respectively.
