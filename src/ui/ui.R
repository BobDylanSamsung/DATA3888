source("src/ui/home.R")
source("src/ui/model_tab.R")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("house")),
    menuItem("EDA", tabName = "eda", icon = icon("magnifying-glass-chart")),
    menuItem("Model", tabName = "model", icon = icon("server")),
    menuItem("Predict", tabName = "predict", icon = icon("chart-line"))
  )
)

create_accordion_box <- function(output_id, title, output_type = c("plot", "table", "text"), height = 500, width = 800, content_after = NULL) {
  output_type <- match.arg(output_type)
  fluidRow(
    box(
      title = title,
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = TRUE,
      width = 12,  # Full-width box
      
      # Output decision based on type
      {
        if (output_type == "plot") {
          plotOutput(output_id, height = height, width = width)
        } else if (output_type == "table") {
          tableOutput(output_id)
        } else {
          verbatimTextOutput(output_id)
        }
      },
      
      # Conditionally render additional content
      if (!is.null(content_after)) {
        content_after
      }
    )
  )
}


eda_tab <- tabItem(tabName = "eda",
  create_accordion_box("expression_histogram", "Expression Histogram", "plot"),
  create_accordion_box("summary_table", "Summary Statistics Table", "table"),
  create_accordion_box("summary_statistics_boxplot", "Summary Statistics Boxplot", "plot"),
  create_accordion_box("expression_per_sample_boxplot", "Expression per Sample Boxplot", "plot"),
  create_accordion_box("pca_plot", "PCA Plot", "plot"),
  create_accordion_box("variable_genes_heatmap", "Heatmap of Top Genes", "plot"),
  plotOutput("variable_genes_heatmap")
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
         column(width = 3, valueBoxOutput("risk_box", width = NULL)),
         column(width = 3, valueBoxOutput("group_box", width = NULL)),
         column(width = 3, valueBoxOutput("hr_box", width = NULL)),
         column(width = 3, valueBoxOutput("percentile_box", width = NULL))
       )
     )
   ),
   
   fluidRow(
     # Top genes box
     box(
       title = "Top Contributing Genes",
       status = "info",
       solidHeader = TRUE,
       width = 6,
       h4("Genes Most Affecting Prediction"),
       tableOutput("top_genes_table"),
       HTML("<p><i>Positive contribution values indicate higher risk, negative values indicate protective effect.</i></p>")
     ),
     # Risk score histogram
     box(
       title = "Risk Score Distribution",
       status = "primary",
       solidHeader = TRUE,
       collapsible = TRUE,
       width = 6,
       plotOutput("pred_risk_hist", height = 350),
       HTML("<p>This histogram shows the distribution of risk scores in the test set. 
  The red line indicates where this patient's risk score falls relative to others.</p>")
     ),
   )
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