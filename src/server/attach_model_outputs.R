attach_model_outputs <- function(output, session) {  # Add session parameter

  # Calculate Brier score and IBS
  brier_results <- calculate_brier_score()
  
  # Render all value boxes
  render_value_boxes(output, brier_results$ibs_value)
  
  # Set up survival analysis data
  km_data <- prepare_survival_data()
  surv_fit <- survfit(Surv(time, status) ~ risk_group, data = km_data)
  surv_diff <- survdiff(Surv(time, status) ~ risk_group, data = km_data)
  
  render_survival_plots(output, km_data, surv_fit, surv_diff)
  
  render_performance_plots(output, brier_results$brier_df)
  render_feature_importance_plot(output)
  render_calibration_plot(output)
}

# Helper function to calculate Brier score
calculate_brier_score <- function() {
  # Define evaluation times for Brier score
  eval_times <- seq(6, floor(max(test_time)), by = 6)
  eval_times <- eval_times[sapply(eval_times, function(t)
    any(test_is_dead == 1 & test_time <= t) && any(test_time >= t))]
  
  # Prepare test data for Score function
  test_score_df <- data.frame(
    time    = test_time,
    is_dead = test_is_dead,
    as.data.frame(X_test_selected_std)
  )
  
  # Calculate Brier score with error handling
  score_res <- tryCatch({
    riskRegression::Score(
      object       = list(CoxBoost = coxboost_model),
      formula      = Surv(time, is_dead) ~ 1,
      data         = test_score_df,
      times        = eval_times,
      metrics      = "brier",
      cens.model   = "km",
      split.method = "none",
      summary      = "ibs"
    )
  }, error = function(e) e)
  
  # Process results
  if (inherits(score_res, "error") || length(eval_times) < 2) {
    warning("Brier-score computation failed: ",
            if (inherits(score_res, "error")) score_res$message else "no eval times")
    brier_df  <- NULL
    ibs_value <- NA_real_
  } else {
    brier_df <- score_res$Brier$score[score_res$Brier$score$model == "CoxBoost", c("times", "Brier")]
    names(brier_df) <- c("time", "brier")
    ibs_value <- score_res$Brier$score[model == "CoxBoost" & times == max(times), IBS]
  }
  
  return(list(brier_df = brier_df, ibs_value = ibs_value, eval_times = eval_times))
}

# Helper function to render value boxes
render_value_boxes <- function(output, ibs_value) {
  output$opt_steps_box <- renderValueBox({
    shinydashboard::valueBox(optimal_steps, "Optimal Boosting Steps",
                             icon = icon("forward-fast"), color = "purple")
  })
  
  output$c_index_box <- renderValueBox({
    shinydashboard::valueBox(round(c_index, 3), "Test-set C-index",
                             icon = icon("chart-line"), color = "green")
  })
  
  output$n_features_box <- renderValueBox({
    shinydashboard::valueBox(length(selected_features), "Features in Final Model",
                             icon = icon("dna"), color = "blue")
  })
  
  output$ibs_box <- renderValueBox({
    # Format IBS for display
    disp <- if (is.null(ibs_value) || length(ibs_value) == 0 || !is.finite(ibs_value)) {
      "—"  # em-dash when not available
    } else {
      format(round(ibs_value, 3), nsmall = 3)
    }
    
    shinydashboard::valueBox(
      value = disp,
      subtitle = "Integrated Brier Score",
      icon = icon("bullseye"),
      color = "orange"
    )
  })
}

# Prepare data for survival analysis
prepare_survival_data <- function() {
  data.frame(
    time = test_time,
    status = test_is_dead,
    risk_group = factor(risk_groups)
  )
}

# Helper function to render survival analysis plots and tables
render_survival_plots <- function(output, km_data, surv_fit, surv_diff) {
  # KM Curve plot
  output$km_curve <- renderPlot({
    p_val <- 1 - pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1)
    p_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
    
    ggsurvplot(
      surv_fit, data = km_data,
      pval = p_text, risk.table = FALSE,
      palette = c("#2E9FDF", "#E7B800"),
      xlab = "Time (months)", ylab = "Survival probability",
      title = "Survival by Predicted Risk Group",
      legend.labs = c("Low Risk", "High Risk"),
      legend.title = "Risk Group",
      ggtheme = theme_bw(), conf.int = TRUE,
      censor.shape = "+", censor.size = 4
    )$plot
  })
  
  # Cumulative hazard curve
  output$chaz_curve <- renderPlot({
    if (length(unique(km_data$risk_group)) < 2)
      return(create_empty_plot("Need ≥2 risk groups"))
    
    ggsurvplot(
      surv_fit, data = km_data, fun = "cumhaz",
      risk.table = FALSE, palette = c("#2E9FDF", "#E7B800"),
      xlab = "Time (months)", ylab = "Cumulative Hazard",
      title = "Cumulative Hazard by Risk Group",
      legend.title = "Risk Group", ggtheme = theme_bw(),
      conf.int = TRUE, censor.shape = "+", censor.size = 4
    )$plot
  })
  
  # Log-rank test results
  output$log_rank_results <- renderTable({
    chisq_value <- surv_diff$chisq
    df <- length(levels(km_data$risk_group)) - 1
    p_value <- 1 - pchisq(chisq_value, df)
    
    data.frame(
      Statistic = c("Chi-Square", "Degrees of Freedom", "P-value", "Interpretation"),
      Value = c(sprintf("%.2f", chisq_value),
                df,
                ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value)),
                ifelse(p_value < 0.05, "Groups differ significantly", "No significant difference"))
    )
  }, rownames = FALSE, striped = TRUE, bordered = TRUE, hover = TRUE)
}

# Helper function to render model performance plots
render_performance_plots <- function(output, brier_df) {
  # ROC curve
  output$roc_curve_plot <- renderPlot({
    timeframe <- 12
    eligible_idx <- which(test_time >= timeframe | test_is_dead == 1)
    
    if (length(eligible_idx) < 10)
      return(create_empty_plot("Insufficient data for ROC"))
    
    binary_outcome <- rep(NA, length(test_is_dead))
    binary_outcome[eligible_idx] <- ifelse(
      test_is_dead[eligible_idx] == 1 & test_time[eligible_idx] <= timeframe, 1, 0
    )
    
    valid_idx <- which(!is.na(binary_outcome))
    
    roc_obj <- pROC::roc(binary_outcome[valid_idx],
                         risk_scores_test[valid_idx],
                         ci = TRUE, quiet = TRUE)
    
    auc_value <- round(as.numeric(roc_obj$auc), 3)
    auc_ci <- round(as.numeric(roc_obj$ci), 3)
    
    pROC::ggroc(roc_obj) +
      geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray") +
      theme_bw() +
      labs(title = paste0(timeframe, "-Month Mortality Prediction"),
           subtitle = paste0("AUC = ", auc_value, " (95% CI: ", auc_ci[1], "-", auc_ci[3], ")"),
           x = "False Positive Rate", y = "True Positive Rate") +
      annotate("text", x = .75, y = .25, label = paste0("n = ", length(valid_idx)), size = 4)
  })
  
  # Time-dependent AUC
  output$time_auc_plot <- renderPlot({
    plot_data <- prepare_time_auc_data()
    if (is.null(plot_data)) 
      return(create_empty_plot("Need events and censored cases or valid times"))
    
    ggplot(plot_data, aes(time, auc)) +
      geom_line(colour = "#2E9FDF", linewidth = 1) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2E9FDF", alpha = .2) +
      geom_hline(yintercept = .5, linetype = "dashed", colour = "grey50") +
      coord_cartesian(ylim = c(0, 1)) +
      scale_x_continuous(breaks = plot_data$eval_times) +
      theme_bw() +
      labs(title = "Time-dependent AUC (test set)",
           x = "Months since diagnosis", y = "AUC")
  })
  
  # Brier curve
  output$brier_curve <- renderPlot({
    if (is.null(brier_df))
      return(create_empty_plot("Not enough follow-up"))
    
    ggplot(brier_df, aes(time, brier)) +
      geom_line(colour = "#2E9FDF", linewidth = 1) +
      geom_hline(yintercept = .25, linetype = "dashed", colour = "grey60") +
      coord_cartesian(ylim = c(0, 1)) +
      scale_x_continuous(breaks = brier_df$time) +
      theme_bw() +
      labs(title = "Time-dependent Brier Score (test set)",
           x = "Months since diagnosis", y = "Brier score")
  })
}

# Prepare data for time-dependent AUC plot
prepare_time_auc_data <- function() {
  delta <- ifelse(test_is_dead == 1, 1, 0)
  T <- test_time
  marker <- risk_scores_test
  
  n_event <- sum(delta == 1)
  n_cens <- sum(delta == 0)
  if (n_event == 0 || n_cens == 0) return(NULL)
  
  eval_times2 <- seq(6, floor(max(T)), by = 6)
  eval_times2 <- eval_times2[sapply(eval_times2, function(t)
    any(delta == 1 & T <= t) && any(delta == 0 & T >= t))]
  if (length(eval_times2) < 2) return(NULL)
  
  td <- timeROC::timeROC(T, delta, as.numeric(marker), cause = 1,
                         times = eval_times2, iid = TRUE,
                         weighting = "marginal")
  
  auc_df <- data.frame(
    time = td$times,
    auc = td$AUC,
    se = td$inference$vect_sd_1
  )
  auc_df$lower <- pmax(auc_df$auc - 1.96 * auc_df$se, 0)
  auc_df$upper <- pmin(auc_df$auc + 1.96 * auc_df$se, 1)
  auc_df$eval_times <- eval_times2
  
  return(auc_df)
}

# Helper to create empty placeholder plot
create_empty_plot <- function(message) {
  ggplot() + 
    annotate("text", x = .5, y = .5, label = message) +
    theme_void()
}

# Add feature importance plot function
render_feature_importance_plot <- function(output) {
  output$feature_importance_plot <- renderPlot({
    # Extract coefficients from the CoxBoost model
    coefs <- coef(coxboost_model)[selected_features]
    
    # Convert to data frame for plotting
    coef_df <- data.frame(
      feature = selected_features,
      coefficient = coefs,
      stringsAsFactors = FALSE
    ) %>%
      # Map feature names to gene symbols when possible
      mutate(
        # Look up original gene ID from the safe name
        original_id = lookup$raw[match(feature, lookup$safe)],
        # Get full gene information
        gene_info = GSE28735$featureData$gene_assignment[match(original_id, rownames(GSE28735$featureData))],
        # Extract gene symbol
        gene_symbol = sapply(gene_info, extract_gene_symbol, USE.NAMES = FALSE),
        # Create display label
        display_name = ifelse(!is.na(gene_symbol), 
                              paste0(gene_symbol, " (", feature, ")"), 
                              feature),
        # Set impact color
        impact = ifelse(coefficient > 0, "Risk (Worse Survival)", "Protective (Better Survival)"),
        # Sort by absolute value for plotting
        abs_coef = abs(coefficient)
      ) %>%
      arrange(desc(abs_coef))
    
    # Create horizontal bar plot
    ggplot(coef_df, aes(x = reorder(display_name, abs_coef), y = coefficient, fill = impact)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Protective (Better Survival)" = "#2E9FDF", 
                                   "Risk (Worse Survival)" = "#E7B800")) +
      coord_flip() +
      labs(
        title = "Feature Importance in CoxBoost Survival Model",
        subtitle = "Features sorted by absolute effect size",
        x = "",
        y = "Coefficient (Effect on Survival)",
        fill = "Impact"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
  })
}

render_calibration_plot <- function(output) {
  output$calibration_plot <- renderPlot({
    time_point <- 6  
    
    # Get predicted survival probabilities at the chosen time point
    time_index <- max(which(coxboost_model$bh_time <= time_point), 1)
    
    # Calculate baseline survival at that time
    baseline_survival <- exp(-coxboost_model$bh_hazard[time_index])
    
    # Calculate predicted survival probabilities for test set
    pred_survival <- baseline_survival^exp(risk_scores_test)
    
    # Create a dataframe with predictions and actual outcomes
    calib_data <- data.frame(
      predicted = pred_survival,
      time = test_time,
      status = test_is_dead
    )
    
    # Create risk groups by quantiles 
    num_groups <- 2
    calib_data$risk_group <- cut(calib_data$predicted, 
                                 breaks = quantile(calib_data$predicted, probs = seq(0, 1, length.out = num_groups + 1)),
                                 labels = 1:num_groups, include.lowest = TRUE)
    
    # Calculate observed survival in each group using Kaplan-Meier
    observed_survival <- data.frame()
    predicted_survival <- data.frame()
    
    for (group in levels(calib_data$risk_group)) {
      group_data <- calib_data[calib_data$risk_group == group, ]
      
      # Skip groups with too few samples
      if (nrow(group_data) < 5) next
      
      # Calculate average predicted survival probability for this group
      avg_pred <- mean(group_data$predicted)
      
      # Calculate observed survival at time_point using Kaplan-Meier
      km_fit <- survfit(Surv(time, status) ~ 1, data = group_data)
      
      # Extract survival at time_point
      surv_at_time <- summary(km_fit, times = time_point)
      
      # If we have data for this time point
      if (length(surv_at_time$surv) > 0) {
        observed <- surv_at_time$surv
        obs_lower <- surv_at_time$lower
        obs_upper <- surv_at_time$upper
        
        # Add to results
        observed_survival <- rbind(observed_survival, 
                                   data.frame(
                                     group = as.numeric(group),
                                     observed = observed,
                                     lower = obs_lower,
                                     upper = obs_upper,
                                     n = nrow(group_data)
                                   ))
        
        predicted_survival <- rbind(predicted_survival, 
                                    data.frame(
                                      group = as.numeric(group),
                                      predicted = avg_pred,
                                      n = nrow(group_data)
                                    ))
      }
    }
    
    # Merge observed and predicted data
    if (nrow(observed_survival) > 0 && nrow(predicted_survival) > 0) {
      calib_plot_data <- merge(observed_survival, predicted_survival, by = "group")
      
      # Create the calibration plot
      ggplot(calib_plot_data, aes(x = predicted, y = observed)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(aes(size = n.x), color = "#2E9FDF", alpha = 0.8) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01, color = "#2E9FDF") +
        scale_size_continuous(name = "Sample size", range = c(3, 10)) +
        labs(
          title = paste0(time_point, "-Month Survival Calibration"),
          subtitle = "Model calibration across different risk groups",
          x = "Predicted survival probability", 
          y = "Observed survival probability"
        ) +
        coord_equal() +
        xlim(0, 1) + 
        ylim(0, 1) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "bottom"
        )
    } else {
      # If we don't have data for the calibration plot
      create_empty_plot("Insufficient follow-up data for calibration plot")
    }
  })
}
