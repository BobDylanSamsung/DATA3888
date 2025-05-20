# Function to parse gene expression CSV and prepare for prediction
parse_csv <- function(csv, expected_genes) {
  df <- read.csv(csv, header = TRUE, stringsAsFactors = FALSE)
  
  # Get column names - first is usually gene ID, second is expression
  gene_col <- colnames(df)[1]
  expr_col <- colnames(df)[2]
  
  # Extract gene IDs and expression values
  genes <- df[[gene_col]]
  values <- as.numeric(df[[expr_col]])
  names(values) <- genes
  
  # Create a complete vector with NAs for missing genes
  complete_values <- setNames(rep(NA, length(expected_genes)), expected_genes)
  
  # Fill in values for genes that are close matches
  for (gene in names(values)) {
    # Check if gene ID without 'X' prefix matches any expected gene
    expected_match <- expected_genes[grepl(paste0("^X", gene, "$"), expected_genes)]
    
    if (length(expected_match) > 0) {
      complete_values[expected_match] <- values[gene]
    }
  }
  
  # Transform into a one-row dataframe suitable for prediction
  newdata <- as.data.frame(t(complete_values))
  
  # Check if we have enough matched genes
  matched_count <- sum(!is.na(complete_values))
  if (matched_count == 0) {
    stop("No gene IDs in the CSV match the expected gene IDs. Please check your file format.")
  } else if (matched_count < 0.5 * length(expected_genes)) {
    warning("Only ", matched_count, " out of ", length(expected_genes), 
            " genes were matched. Results may not be accurate.")
  }
  
  return(newdata)
}

# Function to process prediction data and return results
process_prediction <- function(csv_path, model_params) {
  # Extract model parameters
  coxboost_model <- model_params$coxboost_model
  selected_features <- model_params$selected_features
  train_center <- model_params$train_center
  train_scale <- model_params$train_scale
  cutoff <- model_params$cutoff
  risk_groups <- model_params$risk_groups
  risk_scores_test <- model_params$risk_scores_test
  lookup <- model_params$lookup
  
  # Parse the CSV
  newdata <- parse_csv(csv_path, names(train_center))
  
  # Handle missing values for selected features
  for (feature in selected_features) {
    if (is.na(newdata[1, feature])) {
      # For missing values, use the mean from training data
      newdata[1, feature] <- train_center[feature]
    }
  }
  
  # Extract only the selected features
  newdata_selected <- newdata[, selected_features, drop = FALSE]
  
  # Scale using training parameters
  newdata_scaled <- scale(newdata_selected, 
                          center = train_center, 
                          scale = train_scale)
  
  # Make sure we're using optimal steps for prediction
  optimal_steps <- coxboost_model$stepno
  
  # Calculate prediction metrics using CoxBoost model
  risk_score <- as.vector(predict(coxboost_model, 
                                  newdata = newdata_scaled, 
                                  type = "lp",
                                  at.step = optimal_steps))
  
  risk_group <- ifelse(risk_score > cutoff, "High", "Low")
  hr <- exp(risk_score)
  
  # Estimate survival curve for this patient
  lp <- risk_score
  S0t <- data.frame(time = coxboost_model$bh_time, 
                    hazard = coxboost_model$bh_hazard)
  
  # Calculate survival function S(t) = S0(t)^exp(lp)
  surv_curve <- data.frame(
    time = S0t$time,
    survival = exp(-S0t$hazard * exp(lp))
  )
  
  # Find median survival time (when S(t) crosses 0.5)
  median_idx <- which(surv_curve$survival <= 0.5)[1]
  if(!is.na(median_idx)) {
    median_surv <- surv_curve$time[median_idx]
  } else {
    # If survival never goes below 0.5, use the maximum observed time
    median_surv <- max(surv_curve$time)
  }
  
  # Get CoxBoost coefficients from the model
  coefs_matrix <- as.matrix(coxboost_model$coefficients)
  final_coefs <- coefs_matrix[nrow(coefs_matrix), ]
  names(final_coefs) <- coxboost_model$xnames
  
  # Create gene contributions with proper mapping
  gene_contributions <- data.frame(
    Gene = character(length(selected_features)),
    Safe_Name = selected_features,
    Coefficient = numeric(length(selected_features)),
    Expression = as.numeric(newdata_scaled[1, ]),
    Contribution = numeric(length(selected_features)),
    stringsAsFactors = FALSE
  )
  
  # Map each selected feature to its corresponding coefficient
  for (i in seq_along(selected_features)) {
    feature <- selected_features[i]
    
    # Get the original gene ID from lookup
    original_gene <- lookup$raw[lookup$safe == feature]
    gene_contributions$Gene[i] <- original_gene
    
    # Feature names already have X prefix, no need to add again
    if (feature %in% coxboost_model$xnames) {
      # Found an exact match
      gene_contributions$Coefficient[i] <- final_coefs[feature]
      gene_contributions$Contribution[i] <- final_coefs[feature] * gene_contributions$Expression[i]
    }
  }
  
  # Sort by absolute contribution
  gene_contributions <- gene_contributions[order(abs(gene_contributions$Contribution), decreasing = TRUE), ]
  
  # Calculate risk percentile
  risk_percentile <- ecdf(risk_scores_test)(risk_score) * 100
  
  # Return all results
  return(list(
    risk_score = risk_score,
    risk_group = risk_group,
    hr = hr,
    median_surv = median_surv,
    gene_contributions = gene_contributions,
    surv_curve = surv_curve,
    risk_percentile = risk_percentile
  ))
}

# Function to attach prediction outputs to the UI
attach_prediction_outputs <- function(output, input, pred_results, risk_groups, risk_scores_test, c_index) {
  # Display risk summary
  output$risk_summary <- renderUI({
    # Format risk percentile
    percentile_text <- if(pred_results$risk_percentile > 95) {
      "very high"
    } else if(pred_results$risk_percentile > 75) {
      "high"
    } else if(pred_results$risk_percentile > 50) {
      "above average"
    } else if(pred_results$risk_percentile > 25) {
      "below average"
    } else {
      "low"
    }
    
    div(
      h3("Patient Risk Profile"),
      p(paste0("Risk Score: ", round(pred_results$risk_score, 2))),
      p(paste0("Risk Group: ", pred_results$risk_group)),
      p(paste0("Hazard Ratio: ", round(pred_results$hr, 2))),
      p(paste0("Estimated Median Survival: ", round(pred_results$median_surv, 1), " months")),
      p(HTML(paste0("This patient's risk score is in the <strong>", percentile_text, 
                    "</strong> range (", round(pred_results$risk_percentile), "th percentile).")))
    )
  })
  
  # Plot survival curve
  output$survival_curve <- renderPlot({
    # Only render the vertical line if the slider value exists
    req(input$survival_time_slider)
    
    surv_curve <- pred_results$surv_curve
    selected_time <- input$survival_time_slider
    
    # Find the closest time point in our data
    closest_time_idx <- which.min(abs(surv_curve$time - selected_time))
    closest_time <- surv_curve$time[closest_time_idx]
    selected_prob <- surv_curve$survival[closest_time_idx]
    
    ggplot(surv_curve, aes(x = time, y = survival)) +
      geom_line(color = "steelblue", size = 1.2) +
      # Add dashed vertical line at the selected time point
      geom_vline(xintercept = selected_time, linetype = "dashed", color = "red", size = 1) +
      # Add horizontal line at the corresponding probability 
      geom_hline(yintercept = selected_prob, linetype = "dashed", color = "red", size = 1) +
      labs(x = "Time (months)", y = "Survival Probability", 
           title = "Predicted Survival Curve") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title = element_text(face = "bold"))
  })
  
  # Create time slider input
  output$time_slider_ui <- renderUI({
    surv_curve <- pred_results$surv_curve
    max_time <- max(surv_curve$time)
    
    sliderInput(
      "survival_time_slider", 
      "Select time point (months):",
      min = 0,
      max = max_time,
      value = max_time/4,
      step = 1, 
      width = "100%"
    )
  })
  
  # Display survival probability at selected time
  output$survival_probability <- renderText({
    req(input$survival_time_slider)
    surv_curve <- pred_results$surv_curve
    
    # Find the closest time point in our data
    closest_time_idx <- which.min(abs(surv_curve$time - input$survival_time_slider))
    closest_time <- surv_curve$time[closest_time_idx]
    surv_prob <- surv_curve$survival[closest_time_idx]
    
    paste0("At ", round(closest_time, 1), " months, the estimated survival probability is ", 
           round(surv_prob * 100, 1), "%\n")
  })
  
  
  # Top genes table
  output$gene_table <- DT::renderDT({
    contrib_table <- pred_results$gene_contributions[, c("Gene", "Coefficient", "Expression", "Contribution")]
    colnames(contrib_table) <- c("Gene", "Coefficient", "Standardized Expression", "Contribution to Risk")
    
    # Format the table
    datatable(
      contrib_table, 
      options = list(pageLength = 10, dom = 'ftip'),
      rownames = FALSE
    ) %>%
      formatRound(c("Coefficient", "Standardized Expression", "Contribution to Risk"), 3) %>%
      formatStyle(
        'Contribution to Risk',
        background = styleColorBar(contrib_table$`Contribution to Risk`, 'lightblue'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })
  
  # Risk distribution plot
  output$risk_distribution <- renderPlot({
    bins <- 30
    risk_df <- data.frame(
      score = risk_scores_test,
      group = risk_groups
    )
    
    ggplot() +
      geom_histogram(data = risk_df, aes(x = score, fill = group), bins = bins, alpha = 0.7) +
      geom_vline(xintercept = pred_results$risk_score, linetype = "dashed", size = 1.2, color = "red") +
      annotate("text", x = pred_results$risk_score, y = max(hist(risk_scores_test, bins, plot = FALSE)$counts) * 0.9, 
               label = "Patient", color = "red", hjust = -0.2) +
      scale_fill_manual(values = c("Low" = "green3", "High" = "firebrick")) +
      labs(x = "Risk Score", y = "Count", 
           title = "Risk Score Distribution with Patient Position",
           fill = "Risk Group") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.position = "bottom")
  })
  
  # Clear any previous error message
  output$error_message <- renderUI({ NULL })
}

# Function to clear prediction outputs
clear_prediction_outputs <- function(output, error_message = NULL) {
  output$risk_summary <- renderUI({ NULL })
  output$survival_curve <- renderPlot({ NULL })
  output$gene_table <- DT::renderDT({ NULL })
  output$risk_distribution <- renderPlot({ NULL })
  output$time_slider_ui <- renderUI({ NULL })
  output$survival_probability <- renderText({ NULL })
  
  if (!is.null(error_message)) {
    output$error_message <- renderUI({
      div(
        class = "alert alert-danger",
        h4("Error in prediction:"),
        p(error_message)
      )
    })
  } else {
    output$error_message <- renderUI({ NULL })
  }
}