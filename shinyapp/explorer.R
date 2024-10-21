# =============================
# Interactive sequence explorer
# =============================

# Load required packages
library(shiny)
library(DT)
library(dplyr)

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
  theme = shinytheme("simplex"),
  
  # Allows some JavaScript operations           
  useShinyjs(),

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
    uiOutput("line1_ui"),
    uiOutput("line2_ui"),
    uiOutput("line3_ui")
  )
)

# ==========================================
# Server (instructions for building the app)
# ==========================================

server <- function(input, output) {

  # Filter all possible sequences to the chosen population and lines
  seq <- reactive({
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
  
  output$line1_ui <- renderUI({
    choices <- unique(seq()$V2)
    radioButtons("line1_chosen", "line1", choices)
  })

  output$line2_ui <- renderUI({
    req(input$line1_chosen)
    choices <- unique(seq()[seq()$V2 == input$line1_chosen, "V3"])
    radioButtons("line2_chosen", "line2", choices)
  })

  output$line3_ui <- renderUI({
    req(input$line2_chosen)
    choices <- unique(seq()[seq()$V2 == input$line1_chosen & seq()$V3 == input$line2_chosen, "V4"])
    radioButtons("line3_chosen", "line3", choices)
  })
}

# ====================
# Create the shiny app
# ====================

shinyApp(ui, server)