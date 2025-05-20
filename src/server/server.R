source("src/server/attach_model_outputs.R")
source("src/server/prediction_helpers.R")

attach_eda_outputs <- function(output) {
  # Render Expression Histogram
  output$expression_histogram <- renderPlot({
    hist(gse$eMat, breaks = 100, 
         main = "Distribution of Gene Expression Values", 
         xlab = "Log2 Expression Value", 
         ylab = "Frequency",
         col = "lightblue")
  })
  
  # Render Summary Statistics Table
  output$summary_table <- DT::renderDT({
    DT::datatable(
      summary_per_sample_eda,
      options = list(
        pageLength = 5,
        scrollX = TRUE,
        columnDefs = list(
          # Center align all numeric columns
          list(className = 'dt-center', targets = 1:6),
          # Left align the sample ID column
          list(className = 'dt-left', targets = 0)
        )
      ),
      caption = "Summary statistics for each sample in the dataset",
      rownames = FALSE  # Explicitly disable row names
    ) %>%
      DT::formatRound(columns = c("Min", "Q1", "Median", "Mean", "Q3", "Max"), digits = 2)
  })
  
  # Render Summary Statistics Boxplot
  output$summary_statistics_boxplot <- renderPlot({
    ggplot(summary_long_eda, aes(x = Statistic, y = Value)) +
      geom_boxplot(fill = "skyblue") +
      facet_wrap(~ Statistic, scales = "free") +
      ggtitle("Distribution of Summary Statistics Across Samples") +
      ylab("Expression Value") +
      xlab("") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_blank()
      )
  })
  
  # Render Expression per Sample Boxplot 
  output$expression_per_sample_boxplot <- renderPlot({
    # Create a smaller sample dataset for plotting
    set.seed(42) # For reproducibility
    sample_subset <- sample(colnames(gse$eMat), min(30, ncol(gse$eMat)))
    
    # Count samples by survival status
    alive_count <- sum(gse$phenoData$is_dead == 0)
    deceased_count <- sum(gse$phenoData$is_dead == 1)
    
    sample_data <- data.frame(
      expression = as.vector(gse$eMat[, sample_subset]),
      sample = rep(sample_subset, each = nrow(gse$eMat)),
      status = rep(gse$phenoData$Survival_Status[match(sample_subset, colnames(gse$eMat))], each = nrow(gse$eMat))
    )
    
    ggplot(sample_data, aes(x = sample, y = expression, fill = status)) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
      scale_fill_manual(values = c("Alive" = "skyblue", "Deceased" = "coral")) +
      labs(
        title = "Gene Expression Distribution by Sample",
        subtitle = paste("Showing a subset of samples. Total:", ncol(gse$eMat), 
                         "(Alive:", alive_count, "- Deceased:", deceased_count, ")"),
        y = "Log2 Expression Value",
        x = "Sample ID",
        fill = "Status"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right",
        plot.title = element_text(size = 14)
      )
  })
  
  # Render PCA Plot 
  output$pca_plot <- renderPlot({
    pca_data <- as.data.frame(pca_all$x)
    pca_data$Survival_Status <- gse$phenoData$Survival_Status
    
    # Calculate variance explained
    var_explained <- pca_all$sdev^2 / sum(pca_all$sdev^2) * 100
    
    ggplot(pca_data, aes(x = PC1, y = PC2, color = Survival_Status)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(level = 0.95, linetype = 2) +
      labs(
        title = "Principal Component Analysis of Gene Expression Data",
        subtitle = paste0("PC1 explains ", round(var_explained[1], 1), 
                          "% and PC2 explains ", round(var_explained[2], 1), 
                          "% of variance"),
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
        color = "Survival Status"
      ) +
      scale_color_manual(values = c("Alive" = "blue", "Deceased" = "red")) +
      theme_minimal()
  })
  
  # Render Heatmap of Top Genes 
  output$variable_genes_heatmap <- renderPlot({
    # Set explicit width and height for better rendering
    pheatmap(gse$eMat[top75_genes_eda, ],
             scale = "row",
             show_rownames = FALSE,
             show_colnames = FALSE,
             annotation_col = annotation_col_eda,
             main = "Expression Patterns Across Top 75 Variable Genes",
             annotation_legend = TRUE,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             treeheight_row = 20,
             treeheight_col = 30,
             fontsize = 8,
             cellwidth = NA,
             cellheight = 4)  # Reduced cell height for more genes
  }, height = 400, width = 600)  # Explicit sizing for the plot
}


# Server function for the Shiny app
server <- function(input, output, session) {
  # Initialize reactive values
  pred_made <- reactiveVal(FALSE)
  output$show_model_plots <- reactive({ !pred_made() })
  outputOptions(output, "show_model_plots", suspendWhenHidden = FALSE)
  
  # Attach visualizations to output
  attach_eda_outputs(output)
  attach_model_outputs(output)
  
  # Initialize prediction outputs
  clear_prediction_outputs(output)
  
  # Store model parameters for prediction
  model_params <- list(
    coxboost_model = coxboost_model,
    selected_features = selected_features,
    train_center = attr(X_train_selected_std, "scaled:center"),
    train_scale = attr(X_train_selected_std, "scaled:scale"),
    cutoff = median(risk_scores_test),
    risk_groups = risk_groups,
    risk_scores_test = risk_scores_test,
    lookup = lookup
  )
  
  # Handle prediction button click
  observeEvent(input$predict, {
    req(input$gene_csv)
    pred_made(TRUE)
    updateTabsetPanel(session, "main_tabs", selected = "Prediction")
    
    withProgress(message = 'Processing data...', {
      tryCatch({
        # Process the prediction with all parameters in a single object
        pred_results <- process_prediction(
          input$gene_csv$datapath, 
          model_params
        )
        
        # Attach the prediction outputs to the UI
        attach_prediction_outputs(
          output, 
          pred_results, 
          risk_groups,
          risk_scores_test, 
          c_index
        )
        
      }, error = function(e) {
        clear_prediction_outputs(output, e$message)
      })
    })
  })
  
  # Reset to model plots if file re-uploaded
  observeEvent(input$gene_csv, {
    pred_made(FALSE)
    clear_prediction_outputs(output)
    updateTabsetPanel(session, "main_tabs", selected = "Model")
  })
}