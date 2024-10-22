# =============================
# Interactive sequence explorer
# =============================

# Load required packages
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(DT)

# ======
# Inputs
# ======

# Define list of comparators (labels and values)
comparators = c(
  "Nivolumab" = "nivolumab_monotherapy",
  "Cabozantinib plus nivolumab" = "cabozantinib_plus_nivolumab",
  "Nivolumab plus ipilimumab" = "nivolumab_plus_ipilimumab",
  "Lenvatinib plus pembrolizumab" = "lenvatinib_plus_pembrolizumab",
  "Avelumab plus axitinib" = "avelumab_plus_axitinib",
  "Pazopanib" = "pazopanib",
  "Tivozanib" = "tivozanib",
  "Sunitinib" = "sunitinib",
  "Cabozantinib" = "cabozantinib",
  "Lenvatinib plus everolimus" = "lenvatinib_plus_everolimus",
  "Everolimus" = "everolimus",
  "Axitinib" = "axitinib",
  "Sorafenib" = "sorafenib",
  "Best supportive care" = "BSC")

# Define list of populations (labels and values)
populations = c(
  "Population 1" = "pop1",
  "Population 2" = "pop2",
  "Population 3" = "pop3",
  "Population 4" = "pop4")

# Import table with all possible valid sequences
all_seq <- read.csv("data/valid_sequences.csv")

# =============================================
# User interface (layout and appearance of app)
# =============================================

ui <- fluidPage(
  # Website theme
  # theme = shinytheme("simplex"),
  
  titlePanel("Interactive sequence explorer"),
  mainPanel(
    # Chosen population (pop1-pop4)
    selectInput(inputId = "population",
                label = "Patient population",
                choices = populations),
    
    # Number of treatment lines
    selectInput(inputId = "R_maxlines",
                label = "Max lines within the R model",
                choices = c(3, 4),
                selected = 4),
    
    # Lines of treatment
    uiOutput("line1"),
    uiOutput("line2"),
    uiOutput("line3"),
    uiOutput("line4"),
    uiOutput("line5")
  )
)

# ==========================================
# Server (instructions for building the app)
# ==========================================

server <- function(input, output) {
  
  # Reactive expression to filter data table to chosen population and lines
  seq_start <- reactive({
    # Filter to chosen population
    seq <- all_seq[all_seq$V1 == input$population,]
    # If max lines is 3, remove rows with BSC in V6 (as with four max lines,
    # any combination with 4 treatments will then have BSC as the final
    # treatment, but those with less will have an empty final column)
    if (input$R_maxlines == 3) {
      seq <- seq[seq$V6 != "BSC",]
    }
    return(seq)
  })
  
  # Reactive filtering of dataframe based on chosen treatments
  seq_line1 <- reactive({
    seq <- seq_start()
    seq <- seq[seq$V2 == input$l1_chosen,]
    return(seq)
  })
  seq_line2 <- reactive({
    seq <- seq_line1()
    seq <- seq[seq$V3 == input$l2_chosen,]
    return(seq)
  })
  seq_line3 <- reactive({
    seq <- seq_line2()
    seq <- seq[seq$V4 == input$l3_chosen,]
    return(seq)
  })
  seq_line4 <- reactive({
    seq <- seq_line3()
    seq <- seq[seq$V5 == input$l4_chosen,]
    return(seq)
  })
  
  # Reactive display of possible treatments for each line
  output$line1 <- renderUI({
    l1_values <- unique(seq_start()$V2)
    radioGroupButtons(inputId = "l1_chosen",
                      label = "First line treatment",
                      choices = comparators[comparators %in% l1_values],
                      individual = TRUE,
                      checkIcon = list(
                        yes = icon("ok", lib = "glyphicon")))
  })
  output$line2 <- renderUI({
    l2_values <- unique(seq_line1()$V3)
    radioGroupButtons(inputId = "l2_chosen",
                      label = "Second line treatment",
                      choices = comparators[comparators %in% l2_values],
                      individual = TRUE,
                      checkIcon = list(
                        yes = icon("ok", lib = "glyphicon")))
  })
  output$line3 <- renderUI({
    if (nrow(seq_line2()) > 1) {
      l3_values <- unique(seq_line2()$V4)
      radioGroupButtons(inputId = "l3_chosen",
                        label = "Third line treatment",
                        choices = comparators[comparators %in% l3_values],
                        individual = TRUE,
                        checkIcon = list(
                          yes = icon("ok", lib = "glyphicon")))
    }
  })
  
  output$line4 <- renderUI({
    if (nrow(seq_line3()) > 1) {
      l4_values <- unique(seq_line3()$V5)
      radioGroupButtons(inputId = "l4_chosen",
                        label = "Fourth line treatment",
                        choices = comparators[comparators %in% l4_values],
                        individual = TRUE,
                        checkIcon = list(
                          yes = icon("ok", lib = "glyphicon")))
    }
  })

  output$line5 <- renderUI({
    # Different if statement, as the one used above won't work here
    # Instead, we just check if the final column is blank or not
    if (nrow(seq_line4()) >= 1) {
      if (seq_line4()$V6 != "") {
        l5_values <- unique(seq_line4()$V6)
        radioGroupButtons(inputId = "l5_chosen",
                          label = "Fifth line treatment",
                          choices = comparators[comparators %in% l5_values],
                          individual = TRUE,
                          checkIcon = list(
                            yes = icon("ok", lib = "glyphicon")))
      }
    }
  })
}

# ====================
# Create the shiny app
# ====================

shinyApp(ui, server)