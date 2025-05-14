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

# Function to process prediction data and return results
process_prediction <- function(csv_path, model, test_data, risk_scores) {
  expected_genes <- colnames(X_train)
  newdata <- parse_csv(csv_path, expected_genes)
  
  # Calculate prediction metrics
  risk_score <- predict(model, newx = as.matrix(newdata), type = "link")
  risk_group <- ifelse(risk_score > median(risk_scores), "High", "Low")
  hr <- exp(risk_score)
  
  # Estimate median survival based on risk group
  median_surv <- ifelse(
    risk_group == "High",
    median(test_data$months_survived[test_data$predicted_risk == "High"]),
    median(test_data$months_survived[test_data$predicted_risk == "Low"])
  )
  
  # Calculate gene contributions
  gene_contributions <- data.frame(
    Gene = rownames(coef_cox)[nz_idx_cox],
    Coefficient = as.numeric(coef_cox[nz_idx_cox]),
    Expression = as.numeric(newdata[1, rownames(coef_cox)[nz_idx_cox]]),
    Contribution = as.numeric(coef_cox[nz_idx_cox] * newdata[1, rownames(coef_cox)[nz_idx_cox]])
  )
  gene_contributions <- gene_contributions[order(abs(gene_contributions$Contribution), decreasing = TRUE), ]
  # Add personalized risk metrics
  risk_percentile <- ecdf(risk_scores)(risk_score) * 100
  # Return all results
  return(list(
    risk_score = risk_score,
    risk_group = risk_group,
    hr = hr,
    median_surv = median_surv,
    gene_contributions = gene_contributions,
    newdata = newdata,
    risk_percentile = risk_percentile
  ))
}

# Function to attach prediction outputs to the UI
attach_prediction_outputs <- function(output, pred_results, test_pheno, risk_scores_test, c_index) {
  # Basic risk metrics
  output$risk_box <- renderValueBox({
    valueBox(
      round(pred_results$risk_score, 3),
      "Risk Score",
      icon = icon("chart-line"),
      color = ifelse(pred_results$risk_group == "High", "red", "blue")
    )
  })
  
  output$group_box <- renderValueBox({
    valueBox(
      pred_results$risk_group,
      "Risk Group",
      icon = icon("users"),
      color = ifelse(pred_results$risk_group == "High", "red", "blue")
    )
  })
  
  output$hr_box <- renderValueBox({
    valueBox(
      round(pred_results$hr, 3),
      "Hazard Ratio",
      icon = icon("heartbeat"),
      color = ifelse(pred_results$risk_group == "High", "red", "blue")
    )
  })
  
  output$median_box <- renderValueBox({
    valueBox(
      paste(round(pred_results$median_surv, 1), "mo"),
      "Est. Survival",
      icon = icon("calendar"),
      color = ifelse(pred_results$risk_group == "High", "red", "blue")
    )
  })
  
  # Risk percentile
  output$percentile_box <- renderValueBox({
    # Get the risk percentile
    rp <- round(pred_results$risk_percentile)
    
    # Determine color based on risk percentile (use named colors only)
    box_color <- if(rp > 75) "red" else if(rp > 50) "orange" else if(rp > 25) "yellow" else "blue"
    
    valueBox(
      paste0(rp, "%"),
      "Risk Percentile",
      icon = icon("percent"),
      color = box_color
    )
  })
  
  # Top contributing genes table
  output$top_genes_table <- renderTable({
    head(pred_results$gene_contributions, 10)  # Show top 10 contributing genes
  })
  
  # Risk score histogram
  output$pred_risk_hist <- renderPlot({
    hist(risk_scores_test, breaks=20, main="Risk Score Distribution", 
         xlab="Risk Score", col="lightblue")
    abline(v=pred_results$risk_score, col="red", lwd=2)
    text(pred_results$risk_score, 0, "Patient", pos=4, col="red")
  })
  
}
# Function to clear prediction outputs
clear_prediction_outputs <- function(output, error_message = NULL) {
  if (!is.null(error_message)) {
    output$pred_view <- renderPrint({
      cat("Error in prediction:", error_message, "\n\n",
          "Please ensure your CSV file contains the correct gene identifiers and expression values.")
    })
  } else {
    output$pred_view <- renderPrint({ cat("") })
  }
  
  # Clear all valueBox outputs
  output$risk_box <- renderValueBox({
    valueBox("", "", icon = icon("chart-line"), color = "aqua")
  })
  
  output$group_box <- renderValueBox({
    valueBox("", "", icon = icon("users"), color = "aqua")
  })
  
  output$hr_box <- renderValueBox({
    valueBox("", "", icon = icon("heartbeat"), color = "aqua") 
  })
  
  output$median_box <- renderValueBox({
    valueBox("", "", icon = icon("calendar"), color = "aqua")
  })
  
  output$percentile_box <- renderValueBox({
    valueBox("", "", icon = icon("percent"), color = "aqua")
  })
  
  # Clear other outputs
  output$top_genes_table <- renderTable({ data.frame() })
  output$pred_risk_hist <- renderPlot({})
}


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
  
  # Handle prediction button click
  observeEvent(input$predict, {
    req(input$gene_csv)
    pred_made(TRUE)
    updateTabsetPanel(session, "main_tabs", selected = "Prediction")
    
    tryCatch({
      # Process the prediction
      pred_results <- process_prediction(
        input$gene_csv$datapath, 
        final_model_cox, 
        test_pheno, 
        risk_scores_test
      )
      # Attach the prediction outputs to the UI
      attach_prediction_outputs(
        output, 
        pred_results, 
        test_pheno, 
        risk_scores_test, 
        c_index
      )
      
    }, error = function(e) {
      clear_prediction_outputs(output, e$message)
    })
  })
  
  # Reset to model plots if file re-uploaded
  observeEvent(input$gene_csv, {
    pred_made(FALSE)
    updateTabsetPanel(session, "main_tabs", selected = "model_plots")
  })
}