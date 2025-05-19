source("src/server/attach_model_outputs.R")
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