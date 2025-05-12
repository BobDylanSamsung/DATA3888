ui <- fluidPage(
  titlePanel("Survival prediction Shiny App Based on Lasso-Cox"),
  sidebarLayout(
    sidebarPanel(
      h4("Input the genetic table of the new patient"),
      fileInput("gene_csv", "Upload CSV of genes", accept = ".csv"),
      actionButton("predict", "Start Predict"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        id="main_tabs",
        tabPanel("EDA",
                 h4("Expression Histogram"),
                 plotOutput("eda_hist"),
                 
                 h4("Summary Statistics (first six samples)"),
                 tableOutput("eda_summary_table"),
                 
                 h4("Summary Statistics Boxplot (per Sample)"),
                 plotOutput("eda_boxplot"),
                 
                 h4("Expression per Sample Boxplot"),
                 plotOutput("eda_per_sample_boxplot"),
                 
                 h4("PCA: Tumor vs Normal Samples"),
                 plotOutput("eda_pca"),
                 
                 h4("Top 100 Variable Genes Heatmap"),
                 plotOutput("eda_heatmap")
        ),
        tabPanel("Model", 
                 # Lambda plot, explanations, and stats
                 h4("Lambda Cross-Validation Plot"),
                 tags$div(
                   tags$a(
                     "Show explanation",
                     href = "#",
                     `data-toggle` = "collapse",
                     `data-target` = "#lambda-expl",
                     class = "btn btn-link"
                   ),
                   tags$div(
                     id = "lambda-expl",
                     class = "collapse",
                     lambda_expl
                   )
                 ),
                 plotOutput("lambda_plot"),
                 textOutput("lambda_summary"),
                 
                 tags$hr(),
                 
                 # Coefficient plot, explanation, and stats
                 h4("Lasso-Cox Gene Coefficients"),
                 tags$div(
                   tags$a(
                     "Show explanation",
                     href = "#",
                     `data-toggle` = "collapse",
                     `data-target` = "#coef-expl",
                     class = "btn btn-link"
                   ),
                   tags$div(
                     id = "coef-expl",
                     class = "collapse",
                     coef_expl
                   )
                 ),
                 plotOutput("coef_plot"),
                 textOutput("n_selected_genes"),
                 # Kaplan-Meier Survival Curves by Risk Group Plot
                 h4("Kaplan-Meier Survival Curves by Risk Group"),
                 tags$div(
                   tags$a(
                     "Show explanation",
                     href = "#",
                     `data-toggle` = "collapse",
                     `data-target` = "#km-expl",
                     class = "btn btn-link"
                   ),
                   tags$div(
                     id = "km-expl",
                     class = "collapse",
                     km_expl
                   )
                 ),
                 plotOutput("model_km_plot"),
        ),
        tabPanel("Prediction", 
                 conditionalPanel(
                   "output.show_model_plots == false",
                   verbatimTextOutput("risk_text"),
                   verbatimTextOutput("group_text"),
                   verbatimTextOutput("hr_text"),
                   verbatimTextOutput("median_text"),
                   tableOutput("prob_table"),
                   plotOutput("pred_plot"),
                   plotOutput("pred_risk_hist")   
                 )
        )
      ), width = 9
    )
  )
)