eda_tab <- tabItem(
  tabName = "eda",
  
  # -----------------------------------------------------------------
  # 1. Header section with general information
  # -----------------------------------------------------------------
  fluidRow(
    box(
      title = "Exploratory Data Analysis",
      status = "primary",
      solidHeader = TRUE,
      width = 12,
      p("This section presents visualizations to help understand gene expression patterns in the dataset.
        Note: All samples are tumor samples, colored by survival status (blue = alive, red = deceased).")
    )
  ),
  
  # -----------------------------------------------------------------
  # 2. Expression Distribution Section
  # -----------------------------------------------------------------
  fluidRow(
    # Expression Histogram - left side
    column(width = 6,
           box(
             title = "Expression Histogram",
             status = "info",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("expression_histogram", height = "350px")),
             p("This histogram shows the distribution of log2 expression values across all genes and samples.")
           )
    ),
    
    # Sample Boxplots - right side
    column(width = 6,
           box(
             title = "Expression Boxplots by Sample",
             status = "warning",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("expression_per_sample_boxplot", height = "350px")),
             p("Each boxplot shows gene expression distribution per sample. Blue = alive patients, coral = deceased patients.")
           )
    )
  ),
  
  # -----------------------------------------------------------------
  # 3. Summary Statistics Section
  # -----------------------------------------------------------------
  fluidRow(
    # Summary Table 
    column(width = 6,
           box(
             title = "Sample Statistics Table",
             status = "success",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(DT::DTOutput("summary_table")),
             p("Paginated table of summary statistics for each sample.")
           )
    ),
    
    # Statistics Boxplot - right side
    column(width = 6,
           box(
             title = "Statistics Distribution",
             status = "primary",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("summary_statistics_boxplot", height = "350px")),
             p("Distribution of each summary statistic across all samples.")
           )
    )
  ),
  
  # -----------------------------------------------------------------
  # 4. Data Patterns Section
  # -----------------------------------------------------------------
  fluidRow(
    # PCA Plot - left side
    column(width = 6,
           box(
             title = "Principal Component Analysis",
             status = "danger",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("pca_plot", height = "400px")),
             p("PCA plot showing sample clustering. Points are colored by survival status: blue = alive, red = deceased.")
           )
    ),
    
    # Heatmap - right side
    column(width = 6,
           box(
             title = "Top Gene Expression Patterns",
             status = "warning",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("variable_genes_heatmap", height = "400px")),
             p("Heatmap of the 75 most variable genes. Samples are annotated by survival status.")
           )
    )
  ),
  # -----------------------------------------------------------------
  # 5. Survival Analysis Visualizations
  # -----------------------------------------------------------------
  fluidRow(
    # Survival Time Distribution - left side
    column(width = 6,
           box(
             title = "Survival Time Distribution",
             status = "info",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("survival_time_plot", height = "400px")),
             p("Distribution of survival times (in days), stratified by patient outcome (alive vs. deceased).
             Histograms show the count distribution while the curves show density estimation.")
           )
    ),
    
    # Volcano Plot - right side
    column(width = 6,
           box(
             title = "Differential Expression Analysis",
             status = "primary",
             solidHeader = TRUE,
             width = NULL,
             withDefaultSpinner(plotOutput("volcano_plot", height = "400px")),
             p("Volcano plot showing gene expression differences between alive and deceased patients.
             X-axis: log2 fold change, Y-axis: statistical significance.
             Red points indicate statistically significant genes.")
           )
    )
  ),
  # Add information box at the bottom
  fluidRow(
    box(
      title = "About This Analysis",
      status = "success", 
      solidHeader = TRUE,
      width = 12,
      p("This dataset contains tumor samples only, with patient survival information. The visualizations 
        help identify gene expression patterns that might correlate with survival outcomes.")
    )
  )
)