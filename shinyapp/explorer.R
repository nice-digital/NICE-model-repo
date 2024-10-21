# =============================
# Interactive sequence explorer
# =============================

# Load required packages
library(shiny)
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
  "Best supportive care" = "placebo_BSC")

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

    # First line treatment
    uiOutput("line1"),

    DT::DTOutput("sequences")
  )
)

# ==========================================
# Server (instructions for building the app)
# ==========================================

server <- function(input, output) {

  # Reactive expression to produce data table of valid sequences for each population
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

  # Reactive generation of options for first line
  output$line1 <- renderUI({
    l1_values <- unique(seq()$V2)
    radioButtons(inputId = "l1_chosen",
                 label = "First line treatment",
                 choices = comparators[comparators %in% l1_values],
                 inline=TRUE)
  })

  # Reactive expression to filter seq based on first line choice
  seq2 <- reactive({
    # Get current sequences
    seq2 <- seq()
    # Filter to chosen first line treatment
    seq2 <- seq2[seq2$V2 == input$l1_chosen,]
    return(seq2)
  })

  # Render table of sequences
  output$sequences <-  DT::renderDT(
    {seq2()},
    rownames = FALSE)
}

# ====================
# Create the shiny app
# ====================

shinyApp(ui, server)