source("src/ui/home.R")
source("src/ui/model_tab.R")
source("src/ui/eda_tab.R")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("house")),
    menuItem("EDA", tabName = "eda", icon = icon("magnifying-glass-chart")),
    menuItem("Model", tabName = "model", icon = icon("server")),
    menuItem("Predict", tabName = "predict", icon = icon("chart-line"))
  )
)

predict_tab <- tabItem(tabName = "predict",
  fluidRow(
   box(
     title = "Prediction Input",
     status = "primary",
     solidHeader = TRUE,
     width = 12,
     p("To begin, please upload a CSV file containing the gene expression data of a patient."),
     p("The CSV should contain gene identifiers in the first column and expression values in the second column."),
     fileInput("gene_csv", "Upload Gene Expression CSV", accept = ".csv"),
     actionButton("predict", "Generate Prediction", icon = icon("chart-line"), 
                  class = "btn-primary btn-lg")
   )
  ),
  
  # Error message display area
  uiOutput("error_message"),
  
  # Only show results after prediction is made
  conditionalPanel(
   condition = "output.show_model_plots == false",
   
   # First row: Key metrics in value boxes
   fluidRow(
     box(
       title = "Prediction Summary",
       status = "success",
       solidHeader = TRUE,
       width = 12,
       fluidRow(
         column(width = 2, valueBoxOutput("risk_box", width = NULL)),
         column(width = 2, valueBoxOutput("group_box", width = NULL)),
         column(width = 2, valueBoxOutput("hr_box", width = NULL)),
         column(width = 2, valueBoxOutput("median_box", width = NULL)),
         column(width = 4, valueBoxOutput("percentile_box", width = NULL))
       ),
       uiOutput("risk_summary")
     )
   ),
   
   # Second row: Survival curve and risk distribution
   fluidRow(
     # Survival curve
     box(
       title = "Predicted Survival",
       status = "info",
       solidHeader = TRUE,
       width = 6,
       plotOutput("survival_curve", height = 350),
       HTML("<p>This curve shows the estimated survival probability over time. The red dashed lines indicate the median survival time.</p>")
     ),
     # Risk score histogram
     box(
       title = "Risk Score Distribution",
       status = "primary",
       solidHeader = TRUE,
       width = 6,
       plotOutput("risk_distribution", height = 350),
       HTML("<p>This histogram shows how this patient's risk score compares to the test cohort. The red line indicates the patient's position.</p>")
     )
   ),
   
   # Third row: Gene table and gene contributions
   fluidRow(
     # Gene table box
     box(
       title = "Top Contributing Genes",
       status = "warning",
       solidHeader = TRUE,
       width = 6,
       DT::DTOutput("gene_table"),
       HTML("<p><i>Positive contribution values indicate higher risk, negative values indicate protective effect.</i></p>")
     )
   ),
  )
)

ui <- dashboardPage(
  dashboardHeader(title="Survival Analysis for Pancreatic Cancer"),
  sidebar,
  dashboardBody(
    tabItems(
      home_tab,
      eda_tab,
      model_tab,
      predict_tab
    )
  )
)