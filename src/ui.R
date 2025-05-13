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

home_tab <- tabItem(tabName = "home",
  fluidRow(
    box(plotOutput("plot1", height = 250)),
    
    box(
      title = "Controls",
      sliderInput("slider", "Number of observations:", 1, 100, 50)
    )
  )
)


eda_tab <- tabItem(tabName = "eda",
  create_accordion_box("expression_histogram", "Expression Histogram", "plot"),
  create_accordion_box("summary_table", "Summary Statistics Table", "table"),
  create_accordion_box("summary_statistics_boxplot", "Summary Statistics Boxplot", "plot"),
  create_accordion_box("expression_per_sample_boxplot", "Expression per Sample Boxplot", "plot"),
  create_accordion_box("pca_plot", "PCA Plot", "plot"),
  create_accordion_box("variable_genes_heatmap", "Heatmap of Top Genes", "plot"),
  plotOutput("variable_genes_heatmap")
)

model_tab <- tabItem(tabName = "model",
  fluidRow(
    valueBox(round(best_lambda_cox,2), "Best Lambda"),
    valueBox(round(c_index, 2), "C-Index")
  ),
  create_accordion_box("cv_plot_cox", "Cross-Validation Plot for Cox Model", "plot",content_after = lambda_expl),
  create_accordion_box("coef_plot_cox", "Coefficients Bar Plot for Cox Model", "plot", content_after = coef_expl),
  create_accordion_box("survival_curve", "Kaplan-Meier Survival Curve", "plot", content_after = km_expl),
  create_accordion_box("confusion_matrix_cox", "Confusion Matrix for Cox Model", "text")
)

predict_tab <- tabItem(tabName = "predict",
  div("To begin, please upload a csv file containing the gene expression of a patient"),
  fileInput("gene_csv", "Upload CSV of genes", accept = ".csv"),
  actionButton("predict", "Start Predict"),
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