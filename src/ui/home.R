home_tab <- tabItem(tabName = "home",
  # Welcome Section
  fluidRow(
    box(
      title = "Welcome to the Pancreatic Cancer Survival Analysis Tool", 
      status = "primary", 
      solidHeader = TRUE,
      width = 12,
      tags$div(
        class = "text-center",
        tags$h3("Precision Medicine for Pancreatic Cancer Prognosis", 
                style = "color: #3c8dbc; margin-bottom: 20px;"),
        tags$p("This interactive tool uses gene expression data to predict survival outcomes for pancreatic cancer patients. 
  Upload your patient's data to receive personalized risk assessments and prognostic information.",
               style = "font-size: 16px; margin-bottom: 20px;"),
        
      )
    )
  ),
  
  # Disease Background Section
  fluidRow(
    box(
      title = "Understanding Pancreatic Cancer", 
      status = "info", 
      solidHeader = TRUE,
      width = 6,
      tags$div(
        tags$h4("Disease Background", style = "color: #3c8dbc;"),
        tags$p("Pancreatic cancer is one of the most aggressive malignancies with a 5-year survival rate of only 10%. 
  It's typically diagnosed at advanced stages when surgical intervention is no longer effective."),
        
        tags$h4("Clinical Challenge", style = "color: #3c8dbc; margin-top: 15px;"),
        tags$p("A major challenge in treating pancreatic cancer is predicting patient outcomes and tailoring treatment plans
  to individual patients. Traditional clinical markers alone have limited prognostic value."),
        
        tags$h4("Genomic Approach", style = "color: #3c8dbc; margin-top: 15px;"),
        tags$p("Our model leverages high-throughput gene expression data to identify molecular signatures
  associated with survival outcomes, enabling more personalized treatment decisions.")
      )
    ),
    
    box(
      title = "Model Overview", 
      status = "success", 
      solidHeader = TRUE,
      width = 6,
      tags$div(
        tags$h4("Data Source", style = "color: #00a65a;"),
        tags$p("Our model was trained on gene expression data from a cohort of pancreatic cancer patients with annotated
  clinical outcomes, including survival time and disease status."),
        
        tags$h4("Methodology", style = "color: #00a65a; margin-top: 15px;"),
        tags$p("We employed a penalized Cox proportional hazards model (LASSO) to identify gene expression patterns 
  most strongly associated with patient survival while avoiding overfitting."),
        
        tags$h4("Validation", style = "color: #00a65a; margin-top: 15px;"),
        tags$p("The model was validated using rigorous cross-validation techniques and independent test sets, 
  achieving a concordance index (C-index) that significantly outperforms conventional clinical predictors.")
      )
    )
  ),
  
  fluidRow(
    box(
      title = "Data Preparation Methodology", 
      status = "info", 
      solidHeader = TRUE,
      width = 12,
      tags$div(
        tags$h4("Data Sources", style = "color: #00c0ef;"),
        tags$p("This model is built using gene expression data from two publicly available GEO datasets:"),
        tags$ul(
          tags$li(tags$strong("GSE28735:"), "Contains paired pancreatic tumor and adjacent non-tumor tissue samples with survival data."),
          tags$li(tags$strong("GSE62452:"), "Includes pancreatic ductal adenocarcinoma samples with detailed survival information.")
        ),
        
        tags$h4("Preprocessing Pipeline", style = "color: #00c0ef; margin-top: 15px;"),
        tags$p("Our data preparation involved several key steps:"),
        tags$ol(
          tags$li("Identification of common genes between both datasets to ensure consistent feature space."),
          tags$li("Integration of expression matrices while preserving dataset-specific variation."), 
          tags$li("Extraction and standardization of survival metrics (months survived and mortality status)."),
          tags$li("Filtering to include only tumor samples with complete survival information."),
          tags$li("Creation of cross-validation samples for model testing and validation.")
        ),
        
        tags$div(
          class = "alert alert-warning",
          tags$i(class = "fa fa-lightbulb", style = "margin-right: 8px;"),
          "The final dataset includes ", tags$strong("gene expression profiles from multiple pancreatic tumor samples"), 
          " with matched survival outcomes, enabling robust survival prediction modeling."
        )
      )
    )
  ),
  
  fluidRow(
    box(
      title = "References & Resources", 
      status = "primary", 
      solidHeader = TRUE,
      width = 12,
      tags$div(
        class = "row",
        tags$div(
          class = "col-md-4",
          tags$h4("Key Publications", style = "color: #3c8dbc;"),
          tags$ul(
            tags$li(tags$a(
              href = "https://doi.org/10.1111/cts.13563", 
              "Posta M, Győrffy B. (2023). Analysis of a large cohort of pancreatic cancer transcriptomic profiles to reveal the strongest prognostic factors. Clin Transl Sci.",
              target = "_blank"
            )),
            tags$li(tags$a(
              href = "https://doi.org/10.1038/s41598-020-77220-w",
              "Spooner A, et al. (2020). A comparison of machine learning methods for survival analysis of high-dimensional clinical data for dementia prediction. Sci Rep.",
              target = "_blank"
            )),
            tags$li(tags$a(
              href = "https://doi.org/10.1371/journal.pone.0031507",
              "Zhang G, et al. (2012). DPEP1 inhibits tumor cell invasiveness, enhances chemosensitivity and predicts clinical outcome in pancreatic ductal adenocarcinoma. PLoS One.",
              target = "_blank"
            ))
          )
        ),
        tags$div(
          class = "col-md-4",
          tags$h4("Clinical Resources", style = "color: #3c8dbc;"),
          tags$ul(
            tags$li(tags$a(href = "https://www.cancer.org/cancer/pancreatic-cancer.html", "American Cancer Society")),
            tags$li(tags$a(href = "https://www.cancer.gov/types/pancreatic", "National Cancer Institute"))
          )
        ),
        tags$div(
          class = "col-md-4",
          tags$h4("Technical Documentation", style = "color: #3c8dbc;"),
          tags$ul(
            tags$li(tags$a(href = "https://github.com/BobDylanSamsung/DATA3888", "GitHub Repository", target = "_blank")),
            tags$li(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28735", "GSE28735 Dataset", target = "_blank")),
            tags$li(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452", "GSE62452 Dataset", target = "_blank"))
          )
        )
      )
    )
  ),
  
  # How to Use & Team Section
  fluidRow(
    box(
      title = "How to Use This Tool", 
      status = "warning", 
      solidHeader = TRUE,
      width = 6,
      tags$ol(
        tags$li(tags$b("Explore data:"), "Review the exploratory data analysis (EDA) to understand the dataset characteristics."),
        tags$li(tags$b("Understand the model:"), "Examine the model parameters and performance metrics."),
        tags$li(tags$b("Make predictions:"), "Upload gene expression data in CSV format to generate patient-specific predictions."),
        tags$li(tags$b("Interpret results:"), "Review risk scores, survival estimates, and contributing genetic factors.")
      ),
      tags$div(
        class = "alert alert-info",
        tags$i(class = "fa fa-info-circle", style = "margin-right: 8px;"),
        "Required file format: CSV with gene IDs in column 1 and expression values in column 2."
      )
    ),
    
    box(
      title = "Project Team - Biomed 8", 
      status = "danger", 
      solidHeader = TRUE,
      width = 6,
      tags$div(
        class = "row",
        # Dylan
        tags$div(
          class = "col-md-6",
          tags$div(
            class = "team-member",
            style = "border-left: 3px solid #dd4b39; padding-left: 10px; margin-bottom: 15px;",
            tags$h4("Dylan", style = "margin-top: 0; color: #dd4b39;"),
            tags$p("Mathematics and Computer Science")
          )
        ),
        # Howard
        tags$div(
          class = "col-md-6",
          tags$div(
            class = "team-member",
            style = "border-left: 3px solid #dd4b39; padding-left: 10px; margin-bottom: 15px;",
            tags$h4("Howard", style = "margin-top: 0; color: #dd4b39;"),
            tags$p("Data Science")
          )
        ),
        # Shuhan
        tags$div(
          class = "col-md-6",
          tags$div(
            class = "team-member",
            style = "border-left: 3px solid #dd4b39; padding-left: 10px; margin-bottom: 15px;",
            tags$h4("Shuhan", style = "margin-top: 0; color: #dd4b39;"),
            tags$p("Data Science")
          )
        ),
        # Cherry
        tags$div(
          class = "col-md-6",
          tags$div(
            class = "team-member",
            style = "border-left: 3px solid #dd4b39; padding-left: 10px; margin-bottom: 15px;",
            tags$h4("Cherry", style = "margin-top: 0; color: #dd4b39;"),
            tags$p("Finance and Data Science")
          )
        ),
        # Ashley
        tags$div(
          class = "col-md-6",
          tags$div(
            class = "team-member",
            style = "border-left: 3px solid #dd4b39; padding-left: 10px; margin-bottom: 15px;",
            tags$h4("Ashley", style = "margin-top: 0; color: #dd4b39;"),
            tags$p("Data Science")
          )
        ),
        # Zhantao
        tags$div(
          class = "col-md-6",
          tags$div(
            class = "team-member",
            style = "border-left: 3px solid #dd4b39; padding-left: 10px; margin-bottom: 15px;",
            tags$h4("Zhantao", style = "margin-top: 0; color: #dd4b39;"),
            tags$p("Data Science and Financial Mathematics")
          )
        )
      ),
      tags$div(
        class = "text-center",
        tags$p(tags$em("Biomed 8 - DATA3888"), style = "margin-top: 10px;"),
        tags$p(tags$strong("University of Sydney"), style = "margin-top: 5px;")
      )
    )
  )
)