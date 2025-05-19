attach_model_outputs <- function(output) {
  
  # Main function split into logical sub-functions for better organization
  
  # Calculate Brier score and IBS
  brier_results <- calculate_brier_score()
  
  # Render all value boxes
  render_value_boxes(output, brier_results$ibs_value)
  
  # Set up survival analysis data
  km_data <- prepare_survival_data()
  surv_fit <- survfit(Surv(time, status) ~ risk_group, data = km_data)
  surv_diff <- survdiff(Surv(time, status) ~ risk_group, data = km_data)
  
  # Render survival plots and tables
  render_survival_plots(output, km_data, surv_fit, surv_diff)
  
  # Render model performance plots
  render_performance_plots(output, brier_results$brier_df)
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