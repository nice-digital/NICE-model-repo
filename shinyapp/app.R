library(shiny)
library(shinythemes)

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

# Create user interface
ui <- fluidPage(theme = shinytheme("cosmo"),
  
  # App title
  titlePanel(title = span(img(src="exeter_pentag_small.png", height=60, style = "margin-right: 20px;"),
                          "Exeter Oncology Model: Renal Cell Carcinoma edition"),
             windowTitle = "EOM:RCC"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      selectInput("R_maxlines", "Max lines within the R model",
                  c(1, 2, 3, 4), selected = 4),
      selectizeInput("List_comparators", "Comparator list",
                     choices = comparators, selected = comparators,
                     multiple = TRUE,
                     options = list(plugins = list("remove_button")))

    ),
    
    # Main panel for displaying outputs
    mainPanel(
      
    )
  )
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output) {

}

shinyApp(ui, server)