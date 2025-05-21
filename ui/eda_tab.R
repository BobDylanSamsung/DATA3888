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
      p("This section presents visualizations to help understand gene expression patterns in the dataset. Note: All samples are tumor samples, 
        colored by survival status (blue = alive, red = deceased).")
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
             p("This histogram summarizes the distribution of log₂-transformed gene expression values across all samples. The x-axis represents expression levels 
             (log₂ scale), and the y-axis shows the number of gene–sample pairs. The distribution is unimodal and right-skewed, roughly resembling a normal shape, 
             with most values concentrated between 3.5 and 5.5, indicating that most genes have moderate expression levels. The right tail extends beyond 10, 
             suggesting heterogeneity due to highly expressed genes. Fewer genes fall below 2, indicating a smaller proportion of low-expression genes.
             Overall, the distribution is smooth and lacks obvious outliers, suggesting that the data have been well normalized and batch-corrected, 
             providing a solid foundation for downstream differential analysis and model development.
")
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
             p("This plot shows boxplots of log₂-transformed gene expression values for each sample. Blue boxes represent the alive group and orange boxes 
               represent the deceased group. The boxes (interquartile ranges) and whiskers reflect the median trend and variability of expression within each sample. 
               The median values and spread are generally consistent across samples, with no obvious batch effects or unusual fluctuations. Additionally, 
               the expression distributions of alive and deceased patients largely overlap, indicating balanced and high-quality data across groups.")
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
             p("This table lists descriptive statistics of log₂-transformed gene expression values for each of the 103 samples, including minimum (Min), 
             first quartile (Q1), median, mean, third quartile (Q3), and maximum (Max). From the example data, both the median and mean values across samples 
             are concentrated around 4.1, with Q1 around 2.8, Q3 around 5.5, maximum values close to 13, and minimum values near 1. The overall range and 
             interquartile differences are fairly consistent, with no extreme deviations observed in individual samples. This suggests that after preprocessing, 
             batch effects and outliers have been effectively controlled, and expression distributions across samples are well balanced.
")
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
             p("This plot shows boxplots of the minimum (Min), first quartile (Q1), median, mean, third quartile (Q3), and maximum (Max) values of 
               log₂-transformed gene expression across all samples. The medians of these statistics are roughly centered around 1.0 (Min), 2.9 (Q1), 4.15 (Median), 
               4.30 (Mean), 5.50 (Q3), and 12.8 (Max). The interquartile ranges are narrow, and only a few mild outliers are present. 
               This indicates that overall expression levels and data variability are highly consistent across samples.")
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
             p("This PCA scatter plot visualizes dimensionality reduction of 103 samples based on the first two principal components of 
               gene expression data—PC1 explains 21.6% of the variance and PC2 explains 7.3%. Blue points represent the alive group and red points represent the 
               deceased group. The dashed ellipses outline the clustering range of each group. The samples from both survival groups show strong overlap in the 
               PC1–PC2 space, with no clear separation into distinct clusters, indicating that overall gene expression patterns differ only slightly between the
               two outcomes. However, some alive samples are more spread along the positive direction of PC1, suggesting potential prognostic signals along this axis. 
               In total, PC1 and PC2 together explain about 29% of the total variation.")
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
             p("This heatmap displays the expression patterns of the 75 most variable genes across 103 samples, with survival status indicated at the top 
               (blue: alive, red: deceased). Genes are clustered along the rows, revealing two main expression modules: 
               the gene set on the left shows generally higher expression in deceased patients (red blocks), 
               while the gene set on the right is relatively enriched in alive patients (blue blocks). 
               Sample clustering also groups most deceased patients into one cluster and most alive patients into another, though some overlap remains, 
               indicating that while gene expression patterns related to prognosis are present, there is still notable heterogeneity. 
               Overall, this heatmap visually highlights the association between gene expression and patient survival status.")
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
             p("This plot combines histograms and kernel density curves to compare the distribution of follow-up survival times between alive (blue) and 
               deceased (orange) patients. Survival times for the deceased group are mainly concentrated between 0 and 20 months, peaking around 10–15 months, 
               and decline rapidly over time. In contrast, the alive group shows a more dispersed distribution, with a moderate peak around 20–30 months 
               and a long tail extending beyond 60 months, reflecting right-censoring in patients who remained event-free at the end of follow-up.
               Overall, deceased patients have noticeably shorter survival times, while the alive group shows longer and more heterogeneous survival, 
               providing clear evidence for validating survival models and performing risk stratification.")
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
             p("This volcano plot displays gene expression differences between the alive and deceased groups. 
             The x-axis shows the log₂ fold change (with thresholds of >1 or <–1 indicating biological relevance), 
             and the y-axis shows –log₁₀ adjusted P-values (with a threshold of >1.3 corresponding to P<0.05 for statistical significance). 
             These thresholds are marked with dashed lines. Most genes cluster near the center, 
             suggesting no significant expression difference between the two groups. However, some genes meet both the fold change and significance criteria, 
             highlighted in red, indicating they may be biologically and statistically meaningful.
")
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