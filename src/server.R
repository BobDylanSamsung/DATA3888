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
  
  # Render the cross-validation plot for the CoxBoost model
  output$cv_plot_cox <- renderPlot({
    # Plot CoxBoost CV results
    plot(cv_res$cv.res$logplik ~ seq_len(length(cv_res$cv.res$logplik)),
         type = "b", xlab = "Boosting Steps", ylab = "Cross-validated Log Partial Likelihood",
         main = "Cross-Validation for CoxBoost")
    abline(v = optimal_steps, col = "red", lty = 2)
    text(optimal_steps, min(cv_res$cv.res$logplik), 
         paste("Optimal:", optimal_steps), pos = 4, col = "red")
  })
  
  # Render the feature importance plot from Random Forest
  output$coef_plot_cox <- renderPlot({
    # Extract top features from Random Forest
    feature_imp <- data.frame(
      gene_id = selected_features,
      importance = var_imp$importance[selected_features]
    )
    feature_imp <- feature_imp[order(feature_imp$importance, decreasing = TRUE),]
    
    ggplot(feature_imp, aes(x = reorder(gene_id, importance), y = importance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      theme_minimal() +
      labs(title = "Feature Importance from Random Forest", 
           x = "Gene", y = "Minimal Depth Importance")
  })
  
  # Extract CoxBoost coefficients for display
  output$coxboost_coef_plot <- renderPlot({
    # Get coefficients from CoxBoost model
    coxboost_coefs <- as.numeric(coxboost_model$coefficients)
    names(coxboost_coefs) <- colnames(X_train_selected_std)
    
    # Create data frame for plotting
    coef_df <- data.frame(
      gene_id = names(coxboost_coefs),
      coef = coxboost_coefs
    )
    coef_df <- coef_df[order(abs(coef_df$coef), decreasing = TRUE),]
    
    ggplot(coef_df, aes(x = reorder(gene_id, abs(coef)), y = coef, fill = coef > 0)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_minimal() +
      scale_fill_manual(values = c("red","blue"), labels = c("High-risk","Protective")) +
      labs(title = "CoxBoost Model Coefficients", x = "Gene", y = "Coefficient")
  })
  
  # Render survival curves 
  output$survival_curve <- renderPlot({
    # Create risk groups based on median risk score
    test_risk_groups <- factor(risk_groups, levels = c("Low", "High"))
    
    # Create survival fit
    fit_km <- survfit(Surv(test_time, test_is_dead) ~ test_risk_groups)
    
    # Plot with ggsurvplot
    ggsurvplot(fit_km, data = data.frame(time = test_time, status = test_is_dead, 
                                         risk = test_risk_groups),
               pval = TRUE, risk.table = TRUE, conf.int = TRUE, 
               palette = c("blue", "red"), 
               title = "Test Set Survival Curves by Risk Group")
  })
  
  # Render C-index for model performance evaluation
  output$c_index <- renderText({
    paste("Test set C-index:", round(c_index, 2))
  })
  
  # Render confusion matrix for risk prediction
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
process_prediction <- function(csv_path, coxboost_model, selected_features, train_center, train_scale, cutoff) {
  # Parse the CSV
  newdata <- parse_csv(csv_path, colnames(X))
  
  # Extract only the selected features
  newdata_selected <- newdata[, selected_features, drop = FALSE]
  
  # Scale using training parameters
  newdata_scaled <- scale(newdata_selected, 
                          center = train_center, 
                          scale = train_scale)
  
  # Calculate prediction metrics using CoxBoost model
  risk_score <- predict(coxboost_model, newdata = newdata_scaled, type = "lp")
  risk_group <- ifelse(risk_score > cutoff, "High", "Low")
  hr <- exp(risk_score)
  
  # Estimate median survival based on risk group
  median_surv <- ifelse(
    risk_group == "High",
    median(test_time[risk_groups == "High"]),
    median(test_time[risk_groups == "Low"])
  )
  
  # Get CoxBoost coefficients
  coxboost_coefs <- as.numeric(coxboost_model$coefficients)
  names(coxboost_coefs) <- colnames(X_train_selected_std)
  
  # Calculate gene contributions
  gene_contributions <- data.frame(
    Gene = selected_features,
    Coefficient = coxboost_coefs,
    Expression = as.numeric(newdata_scaled[1, ]),
    Contribution = coxboost_coefs * as.numeric(newdata_scaled[1, ])
  )
  gene_contributions <- gene_contributions[order(abs(gene_contributions$Contribution), decreasing = TRUE), ]
  
  # Add personalized risk metrics
  risk_percentile <- ecdf(risk_scores_test)(risk_score) * 100
  
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
  
  # Store model parameters for prediction
  model_params <- list(
    coxboost_model = coxboost_model,
    selected_features = selected_features,
    train_center = attr(X_train_selected_std, "scaled:center"),
    train_scale = attr(X_train_selected_std, "scaled:scale"),
    cutoff = median(risk_scores_test)
  )
  
  # Handle prediction button click
  observeEvent(input$predict, {
    req(input$gene_csv)
    pred_made(TRUE)
    updateTabsetPanel(session, "main_tabs", selected = "Prediction")
    
    tryCatch({
      # Process the prediction
      pred_results <- process_prediction(
        input$gene_csv$datapath, 
        model_params$coxboost_model, 
        model_params$selected_features,
        model_params$train_center,
        model_params$train_scale,
        model_params$cutoff
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
  
  # Reset to model plots if file re-uploaded
  observeEvent(input$gene_csv, {
    pred_made(FALSE)
    updateTabsetPanel(session, "main_tabs", selected = "model_plots")
  })
}