my_theme <- create_theme(
  adminlte_sidebar(
    dark_bg = "#3c8dbc",
    dark_hover_bg = "lightblue",
    dark_color = "#f7f7f7"
  )
)
source("src/ui/home.R")
source("src/ui/model_tab.R")
source("src/ui/eda_tab.R")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem(" Home", tabName = "home", icon = icon("house")),
    menuItem(" EDA", tabName = "eda", icon = icon("magnifying-glass-chart")),
    menuItem(" Model", tabName = "model", icon = icon("server")),
    menuItem(" Predict", tabName = "predict", icon = icon("chart-line"))
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
    # First row: Prediction Summary box
    fluidRow(
      box(
        title = "Prediction Summary",
        status = "success",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        width = 12,
        withDefaultSpinner(uiOutput("risk_summary"))
      )
    ),
    
    # Second row: Survival curve and risk distribution
    fluidRow(
      # Survival curve
      box(
        title = "Predicted Survival",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        width = 6,
        withDefaultSpinner(plotOutput("survival_curve", height = 350)),
        uiOutput("time_slider_ui"),
        htmlOutput("survival_probability"),
         HTML("<p>This plot shows the model-predicted survival probability for a specific patient over time. 
                           The x-axis represents time (in months), and the y-axis shows the predicted survival probability. 
                           The red dashed lines indicate the predicted survival probability at a specific time point. 
                           The slider allows users to adjust the time point interactively and view survival probabilities at different stages. 
                           Overall, the curve smoothly declines from high to low over time.
</p>")
    ),
      # Risk score histogram
      box(
        title = "Risk Score Distribution",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        width = 6,
        withDefaultSpinner(plotOutput("risk_distribution", height = 350)),
        HTML("<p>This bar chart shows the distribution of risk scores among patients in the test set (green for the low-risk group, 
                             red for the high-risk group), with the current patient's risk score marked by a red dashed line. 
                             Most low-risk patients have scores between –1.3 and 0, while high-risk patients are mainly distributed between 0 and 1.4. 
                             This visualization helps quantify an individual patient's risk and place it in the context of the entire cohort, 
                             making it easier to assess the patient's relative risk level.
</p>")
    )
    ),
    
    # Third row: Gene table and gene contributions
    fluidRow(
      # Gene table box
      box(
        title = "Top Contributing Genes",
        status = "warning",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        width = 6,
        withDefaultSpinner(DT::DTOutput("gene_table")),
       HTML("<p>This table calculates the individual contribution of each gene to the patient's risk score by multiplying the model's 
       regression coefficient for each gene with the patient's standardized gene expression value. A positive value indicates the 
       gene increases risk in this patient, while a negative value suggests a protective effect. The larger the absolute value, 
       the greater the gene's impact on prognosis. Since gene expression profiles vary across patients, 
       the ranking and direction of contributions will differ.</p>")
      ),
      box(
        title = "Gene Expression Radar Chart",
        status = "warning",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        width = 6,
        withDefaultSpinner(plotOutput("radar_chart", height = 400)),
        HTML("<p>This radar chart shows how the patient's gene expression levels compare to the population average for the most influential genes. The population average forms a perfect octagon because the expression values are standardized, with mean=0 for the population. Deviations from this octagon indicate where this patient's expression differs from average.</p>")
      )
    ),
  )
)

ui <- dashboardPage(
  dashboardHeader(title="My Dashboard"),
  sidebar,
  dashboardBody(
    use_theme(my_theme),
    tabItems(
      home_tab,
      eda_tab,
      model_tab,
      predict_tab
    )
  )
)
