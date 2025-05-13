attach_eda_outputs <- function(output) {
  # Render Expression Histogram
  output$expression_histogram <- renderPlot({
    hist(gse$eMat, breaks = 100, main = "All Samples Expression Histogram", xlab = "Expression", col = "lightblue")
  })
  
  # Render Summary Statistics Table
  output$summary_table <- renderTable({
    head(summary_per_sample_eda)  # Display first few rows for preview
  })
  
  # Render Summary Statistics Boxplot
  output$summary_statistics_boxplot <- renderPlot({
    ggplot(summary_long_eda, aes(x = Statistic, y = Value)) +
      geom_boxplot(fill = "skyblue") +
      facet_wrap(~ Statistic, scales = "free") +
      ggtitle("Summary Statistics Boxplot (All Samples)") +
      theme_minimal()
  })
  
  # Render Expression per Sample Boxplot
  output$expression_per_sample_boxplot <- renderPlot({
    boxplot(gse$eMat, outline = FALSE, las = 2, main = "Expression per Sample (All)", col = "gray90")
  })
  
  # Render PCA Plot
  output$pca_plot <- renderPlot({
    autoplot(pca_all, data = gse$phenoData, colour = "Tissue") +
      ggtitle("PCA: Tumor vs Normal Samples") +
      theme_minimal()
  })
  
  # Render Heatmap of Top Genes
  # NOTE: heatmap doesnt show in box for some reason
  #       will only render outside box if window is reshaped
  output$variable_genes_heatmap <- renderPlot({
    pheatmap(gse$eMat[top100_genes_eda, ],
             scale = "row",
             show_rownames = FALSE,
             show_colnames = FALSE,
             annotation_col = annotation_col_eda,
             main = "Top 100 Variable Genes - All Samples")
  })
}

attach_model_outputs <- function(output) {
  
  # Render the cross-validation plot for the Cox model
  output$cv_plot_cox <- renderPlot({
    plot(cvfit_cox, main = "CV for Penalized Cox Model")
  })
  
  # Render the coefficients bar plot for the Cox model
  output$coef_plot_cox <- renderPlot({
    ggplot(df_genes_cox, aes(x = reorder(gene_id, coef), y = coef, fill = coef > 0)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_minimal() +
      scale_fill_manual(values = c("red","blue"), labels = c("High-risk","Protective")) +
      labs(title = "Non-zero Coefficients in Penalized Cox Model", x = "Gene", y = "Coefficient")
  })
  
  # Render survival curves if applicable
  output$survival_curve <- renderPlot({
    fit_km <- survfit(test_surv ~ test_pheno$predicted_risk)
    ggsurvplot(fit_km, data = test_pheno, pval = TRUE, risk.table = TRUE, 
               conf.int = TRUE, palette = c("blue", "red"), 
               title = "Test Set Survival Curves by Risk Group")
  })
  
  # Render C-index for model performance evaluation
  output$c_index <- renderText({
    paste("Test set C-index:", round(c_index, 2))
  })
  
  # Render confusion matrix for risk prediction - analogous structure to logistic results
  output$confusion_matrix_cox <- renderText({
    capture.output(cat("\n==== Confusion Matrix ====\n"), 
                   print(cm_caret), collapse = "\n")
  })
}

parse_csv <- function(csv, expected_genes) {
  df <- read.csv(csv, header = TRUE, stringsAsFactors = FALSE)
  
  genes <- df[[1]]
  values <- as.numeric(df[[2]])
  names(values) <- genes
  
  # Create a complete vector with zeroes for missing genes
  complete_values <- setNames(rep(0, length(expected_genes)), expected_genes)
  complete_values[names(values)] <- values
  
  # Transform into a one-row dataframe suitable for prediction
  newdata <- as.data.frame(t(complete_values))
  colnames(newdata) <- expected_genes
  
  return(newdata)
}

server <- function(input, output, session) {
  pred_result <- reactiveValues()
  pred_made <- reactiveVal(FALSE)
  output$show_model_plots <- reactive({ !pred_made() })
  outputOptions(output, "show_model_plots", suspendWhenHidden = FALSE)
  
  # Attach visualisations to output
  attach_eda_outputs(output)
  attach_model_outputs(output)
  
  observeEvent(input$predict, {
    req(input$gene_csv)
    pred_made(TRUE)
    updateTabsetPanel(session, "main_tabs", selected = "Prediction")
    
    # Define expected genes based on your model training
    expected_genes <- colnames(X_train)  # Use genes from training features
    
    newdata <- parse_csv(input$gene_csv$datapath, expected_genes)  # Parse CSV with expected features
    
    # Initialize prediction result container
    pred_result <- list()
    
    # Predict risk scores using Cox model's linear predictors
    pred_result$risk_score <- predict(final_model_cox, newx = as.matrix(newdata), type = "link")
    pred_result$risk_group <- ifelse(pred_result$risk_score > median(risk_scores_test), "High", "Low")
    pred_result$hr <- exp(pred_result$risk_score)
    
    # Optionally fit survival curves based on new predictions
    pred_result$fit_new <- survfit(Surv(time = gse$phenoData$months_survived, 
                                        event = gse$phenoData$is_dead) ~ pred_result$risk_group, 
                                   data = as.data.frame(gse$phenoData))
    
    pred_result$newdata <- newdata
    
    # Render prediction results
    output$pred_view <- renderPrint({
      print(pred_result)
    })
    
    view(pred_result)  # Use view(pred_result) or equivalent for displaying pred results (ensure view exists)
  })
  
  # Reset to model plots if file re-uploaded
  observeEvent(input$gene_csv, {
    pred_made(FALSE)
    updateTabsetPanel(session, "main_tabs", selected = "model_plots") # adjust accordingly
  })
}