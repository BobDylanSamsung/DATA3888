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

# Function to get survival probability at specific time points
get_survival_probabilities <- function(surv_curve, time_points) {
  result <- sapply(time_points, function(t) {
    # Find the closest time point in our data
    closest_idx <- which.min(abs(surv_curve$time - t))
    surv_prob <- surv_curve$survival[closest_idx]
    return(surv_prob)
  })
  
  # Format as percentages
  formatted <- paste0(round(result * 100, 1), "%")
  
  # Name the results
  names(formatted) <- paste0(time_points, " months")
  
  return(formatted)
}

# Function to attach prediction outputs to the UI
attach_prediction_outputs <- function(output, input, pred_results, risk_groups, risk_scores_test, c_index) {
  # Display risk summary
  # Display risk summary (refactored for better visual appeal)
  output$risk_summary <- renderUI({
    # Format risk percentile
    percentile_text <- if(pred_results$risk_percentile > 95) {
      list(text = "very high", color = "#d9534f")  # Red
    } else if(pred_results$risk_percentile > 75) {
      list(text = "high", color = "#f0ad4e")  # Orange
    } else if(pred_results$risk_percentile > 50) {
      list(text = "above average", color = "#5bc0de")  # Light blue
    } else if(pred_results$risk_percentile > 25) {
      list(text = "below average", color = "#5cb85c")  # Green
    } else {
      list(text = "low", color = "#5cb85c")  # Green
    }
    
    # Define time points for survival estimates
    time_points <- c(6, 12, 24, 36)
    
    # Get probabilities for these time points
    probs <- get_survival_probabilities(pred_results$surv_curve, time_points)
    
    div(
      style = "display: flex; flex-wrap: wrap;",
      
      # Left panel - Patient Risk Profile
      div(
        style = "flex: 1; min-width: 300px; padding-right: 20px;",
        h3(tags$i(class = "fa fa-user-md", style = "margin-right: 10px;"), "Patient Risk Profile"),
        
        # Risk metrics section with styling
        div(
          style = "background-color: #f9f9f9; border-radius: 5px; padding: 15px; margin-bottom: 15px;",
          
          div(style = "display: flex; margin-bottom: 12px;",
              div(style = "width: 170px; font-weight: bold;", "Risk Score:"),
              div(style = "flex: 1;", round(pred_results$risk_score, 2))
          ),
          
          div(style = "display: flex; margin-bottom: 12px;",
              div(style = "width: 170px; font-weight: bold;", "Risk Group:"),
              div(style = "flex: 1;", 
                  tags$span(
                    style = paste0("font-weight: bold; color: ", 
                                   ifelse(pred_results$risk_group == "High", "#d9534f", "#5cb85c")),
                    pred_results$risk_group
                  )
              )
          ),
          
          div(style = "display: flex; margin-bottom: 12px;",
              div(style = "width: 170px; font-weight: bold;", "Hazard Ratio:"),
              div(style = "flex: 1;", round(pred_results$hr, 2))
          ),
          
          div(style = "display: flex; margin-bottom: 12px;",
              div(style = "width: 170px; font-weight: bold;", "Estimated Median Survival:"),
              div(style = "flex: 1;", paste0(round(pred_results$median_surv, 1), " months"))
          ),
          
          # Percentile indicator with custom styling
          hr(),
          div(
            style = "margin-top: 10px;",
            p(
              "This patient's risk score is in the ",
              tags$span(
                style = paste0("font-weight: bold; color: ", percentile_text$color),
                percentile_text$text
              ),
              " range (",
              tags$span(style = "font-weight: bold", paste0(round(pred_results$risk_percentile), "th")),
              " percentile)."
            ),
            
            # Add a visual percentile bar
            div(
              style = "height: 8px; background-color: #e9ecef; border-radius: 4px; margin-top: 10px;",
              div(
                style = paste0(
                  "width: ", round(pred_results$risk_percentile), "%; ",
                  "height: 100%; ",
                  "background-color: ", percentile_text$color, "; ",
                  "border-radius: 4px;"
                )
              )
            ),
            div(
              style = "display: flex; justify-content: space-between; font-size: 0.8em; color: #6c757d; margin-top: 3px;",
              div("Lower Risk"),
              div("Higher Risk")
            )
          )
        )
      ),
      
      # Right panel - Survival Probabilities
      div(
        style = "flex: 1; min-width: 300px;",
        h3(tags$i(class = "fa fa-chart-line", style = "margin-right: 10px;"), "Survival Probability Estimates"),
        
        # Create an HTML table with styling
        div(
          style = "background-color: #f9f9f9; border-radius: 5px; padding: 15px;",
          tags$table(
            class = "table table-hover",
            style = "margin-bottom: 0;",
            tags$thead(
              tags$tr(
                tags$th(style = "border-top: none;", "Time Point"),
                tags$th(style = "border-top: none; text-align: right;", "Survival Probability")
              )
            ),
            tags$tbody(
              lapply(seq_along(probs), function(i) {
                prob_value <- as.numeric(sub("%", "", probs[i])) / 100
                
                # Define color based on probability value
                color <- if(prob_value >= 0.7) {
                  "#5cb85c"  # Green for high probability
                } else if(prob_value >= 0.4) {
                  "#f0ad4e"  # Orange for medium
                } else {
                  "#d9534f"  # Red for low probability
                }
                
                tags$tr(
                  tags$td(names(probs)[i]),
                  tags$td(
                    style = "text-align: right;",
                    div(
                      style = "display: flex; align-items: center; justify-content: flex-end;",
                      div(
                        style = paste0(
                          "width: ", round(prob_value * 100), "%; ",
                          "height: 8px; ",
                          "background-color: ", color, "; ",
                          "border-radius: 4px; ",
                          "margin-right: 10px;"
                        )
                      ),
                      tags$span(
                        style = paste0("font-weight: bold; color: ", color),
                        probs[i]
                      )
                    )
                  )
                )
              })
            )
          )
        )
      )
    )
  })
  
  # Add survival probability table
  output$survival_table <- renderUI({
    # Define time points
    time_points <- c(6, 12, 24, 36)
    
    # Get probabilities for these time points
    probs <- get_survival_probabilities(pred_results$surv_curve, time_points)
    
    # Create an HTML table
    tagList(
      h4("Survival Probability Estimates"),
      tags$table(
        class = "table table-bordered table-striped",
        tags$thead(
          tags$tr(
            tags$th("Time Point"),
            tags$th("Survival Probability")
          )
        ),
        tags$tbody(
          lapply(seq_along(probs), function(i) {
            tags$tr(
              tags$td(names(probs)[i]),
              tags$td(probs[i])
            )
          })
        )
      )
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
    # Get gene contributions
    contrib_table <- pred_results$gene_contributions[, c("Gene", "Coefficient", "Expression", "Contribution")]
    
    # Add gene symbols column using similar approach as in feature_importance_plot
    contrib_table$Gene_Symbol <- sapply(contrib_table$Gene, function(gene_id) {
      # Look up original gene ID (it's already in the Gene column, but keeping the code pattern for consistency)
      original_id <- gene_id
      
      # Get full gene information (assumes GSE28735 is available in this scope)
      gene_info <- GSE28735$featureData$gene_assignment[match(original_id, rownames(GSE28735$featureData))]
      
      # Extract gene symbol
      gene_symbol <- extract_gene_symbol(gene_info)
      
      # Return NA or placeholder if no symbol found
      if(is.na(gene_symbol) || gene_symbol == "") {
        return(NA)
      } else {
        return(gene_symbol)
      }
    })
    
    # Reorder columns to put Gene_Symbol first, followed by Gene (now Gene ID)
    contrib_table <- contrib_table[, c("Gene_Symbol", "Gene", "Coefficient", "Expression", "Contribution")]
    
    # Rename columns
    colnames(contrib_table) <- c("Gene Symbol", "Gene ID", "Coefficient", "Standardized Expression", "Contribution to Risk")
    
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
      ) %>%
      formatStyle(
        'Gene Symbol',
        fontWeight = 'bold'  # Make gene symbols bold for emphasis
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
  
  # Add radar chart for gene expression comparison
  output$radar_chart <- renderPlot({
    # Get all genes with non-zero coefficients
    nonzero_coef_genes <- pred_results$gene_contributions %>% 
      filter(Coefficient != 0) %>%
      pull(Safe_Name)
    
    # Select the genes to display - prioritize non-zero coefficient genes
    if(length(nonzero_coef_genes) >= 8) {
      # If we have 8+ nonzero genes, use those
      genes_to_use <- nonzero_coef_genes[1:8]
    } else {
      # Otherwise use all nonzero genes plus top contributing genes to reach 8
      remaining_genes <- setdiff(pred_results$gene_contributions$Safe_Name, nonzero_coef_genes)
      remaining_needed <- 8 - length(nonzero_coef_genes)
      genes_to_use <- c(nonzero_coef_genes, head(remaining_genes, remaining_needed))
    }
    
    # Get the data for selected genes
    selected_genes_data <- pred_results$gene_contributions %>%
      filter(Safe_Name %in% genes_to_use)
    
    # Extract gene IDs and standardized expression values
    genes <- selected_genes_data$Safe_Name
    patient_expr <- selected_genes_data$Expression
    
    # Population average is 0 for standardized data
    pop_avg <- rep(0, length(genes))
    
    # Get min and max for scaling
    min_expr <- min(c(patient_expr, pop_avg)) - 0.5
    max_expr <- max(c(patient_expr, pop_avg)) + 0.5
    
    # Get gene symbols for labels
    gene_symbols <- sapply(genes, function(feature) {
      # Look up original gene ID from the safe name
      original_id <- lookup$raw[match(feature, lookup$safe)]
      
      # Get full gene information
      gene_info <- GSE28735$featureData$gene_assignment[match(original_id, rownames(GSE28735$featureData))]
      
      # Extract gene symbol
      gene_symbol <- extract_gene_symbol(gene_info)
      
      # Return gene symbol or feature name if symbol not available
      ifelse(!is.na(gene_symbol) && gene_symbol != "", gene_symbol, feature)
    })
    
    # Create dataframe for radar chart
    radar_data <- data.frame(
      # First row is max values, second row is min values
      rbind(
        rep(max_expr, length(genes)),  # Max values
        rep(min_expr, length(genes)),  # Min values
        pop_avg,                       # Population average (0 for standardized data)
        patient_expr                   # Patient values
      )
    )
    
    # Set column names using gene symbols
    colnames(radar_data) <- gene_symbols
    
    # Set row names for clarity
    rownames(radar_data) <- c("Max", "Min", "Population_Average", "Patient")
    
    # Create radar chart using fmsb package
    library(fmsb)
    
    # Set up the plot
    par(mar = c(2, 2, 3, 2))
    
    # Create radar chart
    radarchart(
      radar_data,
      pfcol = c(NA, rgb(0, 0.6, 0.8, 0.3)),     # Population average area color
      pcol = c(1, 1, rgb(0, 0.6, 0.8), "red"),  # Line colors
      plwd = c(1, 1, 3, 3),                      # Line widths
      plty = c(1, 1, 1, 1),                      # Line types
      cglcol = "grey",                           # Grid color
      cglty = 1,                                 # Grid line type
      axislabcol = "grey30",                     # Axis label color
      caxislabels = seq(round(min_expr, 1), round(max_expr, 1), length.out = 5),  # Custom axis labels
      cglwd = 0.8,                               # Grid line width
      title = "Gene Expression Relative to Population"
    )
    
    # Add a legend
    legend(
      "topright",
      legend = c("Patient", "Population Average"),
      col = c("red", rgb(0, 0.6, 0.8)),
      lwd = 3,
      cex = 0.8,
      bty = "n"
    )
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
  output$radar_chart <- renderPlot({ NULL })
  output$survival_table <- renderUI({ NULL })
  
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