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
     ),
     # Gene contribution visualization
     box(
       title = "Gene Contribution Visualization",
       status = "success",
       solidHeader = TRUE,
       width = 6,
       plotOutput("gene_contributions", height = 400),
       HTML("<p>This plot shows the top genes contributing to the risk prediction. Red bars increase risk while blue bars decrease risk.</p>")
     )
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