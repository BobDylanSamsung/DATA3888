## Cox Model Shiny App (app.R)
# Load required libraries
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggfortify)
library(knitr)
library(tidyverse)
library(dplyr)
library(survival)
library(glmnet)
library(survminer)
library(caret)
library(shiny)
library(GEOquery)

par(mar = c(4, 3, 3, 1) + 0.1)

##################################################
## Text and outputs
##################################################
lambda_expl <- HTML('
  <b>Horizontal axis:</b> Log(λ). When λ gets larger (to the left), the penalty becomes stronger and the model becomes more sparse. When λ gets smaller (to the right), the penalty becomes weaker and more genes are kept.<br>
  <b>Vertical axis:</b> Partial Likelihood Deviance. A smaller value means the model fits better.<br>
  Each red dot shows the average deviance for a λ value. The gray lines show the range of error across folds.<br>
  The left dashed line shows λ.min, where the error is the smallest. The right dashed line shows λ.1se.<br>
  We chose λ.min for our model (see below).<br>
')

coef_expl <- HTML('
  The absolute value of each gene\'s coefficient reflects its contribution to the risk score: the larger the value, the more important for distinguishing prognosis.<br>
  <span style="color:blue"><b>Blue bars</b></span> indicate protective genes (negative coefficient): higher expression, lower risk.<br>
  <span style="color:red"><b>Red bars</b></span> indicate high-risk genes (positive coefficient): higher expression, higher risk.<br>
')

km_expl <- HTML('
  Risk grouping: Patients were divided into "low-risk group" (red) and "high-risk group" (cyan) based on the median risk score of all samples. <br>
  Kaplan–Meier curve: The horizontal axis shows follow-up time (months), and the vertical axis shows survival probability.<br>
  The red curve stays above the cyan curve at all times, indicating that patients in the low-risk group have a much higher survival rate compared to those in the high-risk group.<br>
  Log-rank test p-value: p < 0.0001, showing that the difference between the two survival curves is statistically highly significant. This means that the risk score based on the 29 selected genes can clearly distinguish between good and poor prognosis.<br>
  Risk table (Number at risk): Shown below the curves, it displays the number of patients still under follow-up at each time point, helping to evaluate the reliability of survival estimates in later stages.
')

# ----------------------------------------------------------
# Output rendering functions for use in server()
# ----------------------------------------------------------

# Lambda cross-validation plot
lambda_plot_output <- function(cvfit) {
  renderPlot({
    plot(cvfit, main = "Cross-validation for Penalized Cox Model")
  })
}

# Lambda summary text
lambda_summary_output <- function(best_lambda, gene_coefs) {
  renderText({
    paste0(
      "Best λ by cross-validation (λ.min): ", signif(best_lambda, 5), ". ",
      "This value selected ", nrow(gene_coefs), " predictors (genes) for the final model."
    )
  })
}

# Lasso-Cox gene coefficients barplot
coef_plot_output <- function(gene_coefs) {
  renderPlot({
    ggplot(gene_coefs, aes(x = reorder(gene, coef), y = coef, fill = coef > 0)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c("red", "blue"),
                        labels = c("high-risk gene", "protective gene")) +
      labs(x = NULL, y = "Coefficient", title = "Lasso-Cox gene coefficients") +
      theme_minimal()
  })
}

# Number of selected gene summary
n_selected_genes_output <- function(gene_coefs) {
  renderText({
    ng <- nrow(gene_coefs)
    paste0("After training, a total of ", ng, " genes were selected. ",
           "The model has ", sum(gene_coefs$coef > 0), " high-risk (positive) and ",
           sum(gene_coefs$coef < 0), " protective (negative) genes."
    )
  })
}

# Kaplan-Meier survival plot for model
model_km_plot_output <- function(fit_km, train_df) {
  renderPlot({
    survminer::ggsurvplot(fit_km, data = train_df, pval = TRUE, risk.table = TRUE,
                          palette = c("red", "cyan"),
                          title = "Kaplan-Meier Survival Curves by Risk Group")
  })
}

# Single patient predicted KM plot
pred_plot_output <- function(pred_result) {
  renderPlot({
    req(pred_result$fit_new)
    plot(pred_result$fit_new,
         xlab = "Months",
         ylab = "Predicted survival probability",
         main = "Patient's Predicted Survival Curve",
         col = "blue", lwd = 2)
    abline(h = 0.5, col = "red", lty = 2)
  })
}

# Patient risk score on cohort histogram
pred_risk_hist_output <- function(train_lp, pred_result) {
  renderPlot({
    req(pred_result$risk_score)
    hist(train_lp, breaks = 30, xlab = "Risk score",
         main = "Risk score distribution (cohort)", col = "lightgray")
    abline(v = as.numeric(pred_result$risk_score), col = "red", lwd = 3)
    legend("topright", legend = "This patient", col = "red", lwd = 2)
  })
}

# Patient-specific outputs
risk_text_output <- function(pred_result) {
  renderText({
    req(pred_result$risk_score)
    paste0("Risk Score: ", round(pred_result$risk_score, 3))
  })
}

group_text_output <- function(pred_result) {
  renderText({
    req(pred_result$risk_group)
    paste0("Risk Group (High/Low): ", pred_result$risk_group)
  })
}

hr_text_output <- function(pred_result) {
  renderText({
    req(pred_result$hr)
    paste0("Hazard ratio: ", round(pred_result$hr, 3))
  })
}

median_text_output <- function(pred_result) {
  renderText({
    req(pred_result$fit_new)
    idx <- which.min(abs(pred_result$fit_new$surv - 0.5))
    med_time <- round(pred_result$fit_new$time[idx], 1)
    paste0("Predicted median survival time: ", med_time, " months")
  })
}

prob_table_output <- function(pred_result) {
  renderTable({
    req(pred_result$fit_new)
    times <- c(6, 12, 24, 36)
    surv_summary <- summary(pred_result$fit_new, times = times)
    data.frame(
      Month = surv_summary$time,
      Survival_Prob = round(surv_summary$surv, 3)
    )
  })
}

eda_hist_output <- function(gse_eda) {
  renderPlot({
    hist(as.vector(gse_eda$eMat), breaks = 100, 
         main = "All Samples Expression Histogram", 
         xlab = "Expression", col = "lightblue")
  })
}

eda_summary_table_output <- function(summary_per_sample_eda) {
  renderTable({
    head(summary_per_sample_eda, 6)
  }, rownames = TRUE)
}

eda_boxplot_output <- function(summary_long_eda) {
  renderPlot({
    ggplot(summary_long_eda, aes(x = Statistic, y = Value)) +
      geom_boxplot(fill = "skyblue") +
      facet_wrap(~ Statistic, scales = "free") +
      ggtitle("Summary Statistics Boxplot (All Samples)") +
      theme_minimal()
  })
}

eda_per_sample_boxplot_output <- function(gse_eda) {
  renderPlot({
    boxplot(gse_eda$eMat, outline = FALSE, las = 2, 
            main = "Expression per Sample (All)", col = "gray90")
  })
}

eda_pca_output <- function(pca_all, gse_eda) {
  renderPlot({
    autoplot(pca_all, data = gse_eda$phenoData, colour = "Tissue") +
      ggtitle("PCA: Tumor vs Normal Samples") +
      theme_minimal()
  })
}

eda_heatmap_output <- function(gse_eda, top100_genes_eda, annotation_col_eda) {
  renderPlot({
    pheatmap(gse_eda$eMat[top100_genes_eda, ],
             scale = "row",
             show_rownames = FALSE,
             show_colnames = FALSE,
             annotation_col = annotation_col_eda,
             main = "Top 100 Variable Genes - All Samples")
  })
}

# Function to attach all output renderings in one call
attach_render_outputs <- function(output, 
         cvfit, 
         best_lambda, 
         gene_coefs, 
         fit_km, 
         train_df, 
         train_lp, 
         pred_result,
         gse_eda,
         summary_per_sample_eda,
         summary_long_eda,
         pca_all,
         top100_genes_eda,
         annotation_col_eda) {
  
  # Model tab outputs
  output$lambda_plot      <- lambda_plot_output(cvfit)
  output$lambda_summary   <- lambda_summary_output(best_lambda, gene_coefs)
  output$coef_plot        <- coef_plot_output(gene_coefs)
  output$n_selected_genes <- n_selected_genes_output(gene_coefs)
  output$model_km_plot    <- model_km_plot_output(fit_km, train_df)
  
  # Prediction tab outputs
  output$pred_plot        <- pred_plot_output(pred_result)
  output$pred_risk_hist   <- pred_risk_hist_output(train_lp, pred_result)
  output$risk_text        <- risk_text_output(pred_result)
  output$group_text       <- group_text_output(pred_result)
  output$hr_text          <- hr_text_output(pred_result)
  output$median_text      <- median_text_output(pred_result)
  output$prob_table       <- prob_table_output(pred_result)
  
  # EDA tab outputs
  output$eda_hist               <- eda_hist_output(gse_eda)
  output$eda_summary_table      <- eda_summary_table_output(summary_per_sample_eda)
  output$eda_boxplot            <- eda_boxplot_output(summary_long_eda)
  output$eda_per_sample_boxplot <- eda_per_sample_boxplot_output(gse_eda)
  output$eda_pca                <- eda_pca_output(pca_all, gse_eda)
  output$eda_heatmap            <- eda_heatmap_output(gse_eda, top100_genes_eda, annotation_col_eda)
  
  invisible(output)
}
