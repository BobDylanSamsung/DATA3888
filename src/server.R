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

attach_model_ouputs <- function(output) {
  # Render the cross-validation plot
  output$cv_plot <- renderPlot({
    plot(cvfit_binom, main = "CV for Lasso-Logistic Model")
  })
  
  # Render the coefficients bar plot
  # NOTE: plot doesnt show in box for some reason
  #       will only render outside box if window is reshaped
  output$coef_plot <- renderPlot({
    ggplot(df_genes, aes(x = reorder(gene_id, coef), y = coef)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      coord_flip() +
      theme_minimal() +
      labs(title = "Non-zero Coefficients in Lasso-Logistic Model")
  })
  
  # Render confusion matrix
  output$confusion_matrix <- renderText({
    paste(capture.output(cm_results), collapse = "\n")
  })
}

parse_csv <- function(csv) {
  df <- read.csv(csv, header = TRUE)
  
  genes <- df[[1]]
  values <- as.numeric(df[[2]])
  names(values) <- genes
  
  newdata <- as.data.frame(t(values))
  colnames(newdata) <- genes
  
  return(newdata)
}

server <- function(input, output, session) {
  pred_result <- reactiveValues()
  # Flag for whether we are showing the model or prediction
  pred_made <- reactiveVal(FALSE)
  output$show_model_plots <- reactive({ !pred_made() })
  outputOptions(output, "show_model_plots", suspendWhenHidden = FALSE)
  
  # attach visualisations to output
  attach_eda_outputs(output)
  attach_model_ouputs(output)
  
  observeEvent(input$predict, {
    req(input$gene_csv)
    pred_made(TRUE)
    updateTabsetPanel(session, "main_tabs", selected = "Prediction")
    
    newdata <- parse_csv(input$gene_csv$datapath)
    
    cat(dim(newdata))
    
    # Assuming newdata has been correctly prepared (as a one-row dataframe with gene IDs in columns)
    pred_result <- list()  # Initialize a list to store prediction results
    # Predict the class probabilities for the newdata
    pred_probabilities <- predict(final_model_binom, newx = as.matrix(newdata), type = "response")
    cutoff <- 0.5
    pred_result$risk_group <- ifelse(pred_probabilities > cutoff, "High", "Low")
    pred_result$probabilities <- pred_probabilities
    view(pred_result)
  })
  
  # Reset to model plots if file re-uploaded
  observeEvent(input$gene_csv, pred_made(FALSE))
}