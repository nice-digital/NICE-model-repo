library(shiny)
library(shinythemes)
library(data.table)
library(gtools)

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

# Import the sequencing functions
f_path = "../3_Functions"
source(file.path(f_path, "sequencing/sequences.R"))

# Create user interface (layout and appearance of app)
ui <- fluidPage(theme = shinytheme("cosmo"),
                
  # App title
  titlePanel(title = span(img(src="exeter_pentag_small.png", height=60, style = "margin-right: 20px;"),
                          "Exeter Oncology Model: Renal Cell Carcinoma edition"),
             windowTitle = "EOM:RCC"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      #selectInput("i_nr_population", "Number of populations",
      #            c(1, 2, 3, 4), selected = 4),
      selectInput("R_maxlines", "Max lines within the R model",
                  c(1, 2, 3, 4), selected = 4),
      selectizeInput("List_comparators", "Comparator list",
                     choices = comparators, selected = c("nivolumab_monotherapy",
                                                         "cabozantinib_plus_nivolumab",
                                                         "nivolumab_plus_ipilimumab",
                                                         "lenvatinib_plus_pembrolizumab",
                                                         "avelumab_plus_axitinib"),
                     multiple = TRUE,
                     options = list(plugins = list("remove_button"))),
      actionButton("seq_button", "Find possible treatment sequences")

    ),
    
    # Main panel for displaying outputs
    mainPanel(
      textOutput("selected_var"),
      tableOutput("sequences")
    )
  )
)

# Instructions for building the app
server <- function(input, output) {

  # Line of text describing number of lines and treatment
  output$selected_var <- renderText({
    paste0("Model with ", input$R_maxlines, " lines of treatment and ",
           length(input$List_comparators), " treatment options")}) |>
    bindEvent(input$seq_button)

  # Table of all possible sequences
  # Only render when click button
  output$sequences <- renderTable(head(f_generate_sequences(
    comparators = input$List_comparators, maxlines = as.numeric(input$R_maxlines))), striped=TRUE) |>
    bindEvent(input$seq_button)

}

# Create the shiny app
shinyApp(ui, server)