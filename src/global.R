# Load required libraries
library(caret)
library(CoxBoost)
library(dplyr)
library(DT)
library(future)
library(fmsb)
library(GEOquery)
library(ggplot2)
library(ggfortify)
library(glmnet)
library(ggrepel)
library(Hmisc)
library(knitr)
library(limma)
library(pheatmap)
library(promises)
library(randomForestSRC)
library(pec)
library(pROC)
library(prodlim)
library(RColorBrewer)
library(riskRegression)
library(rms)
library(reshape2)
library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(survival)
library(survminer)
library(tidyverse)
library(timeROC)
library(fresh)

set.seed(3888)
plan(multisession)

## ================================================================
##  S3 method required by riskRegression::Score 
## ================================================================
predictRisk.CoxBoost <- function(object, newdata, times, ...) {
  
  # make sure the covariate order is identical
  X <- as.matrix(newdata[, object$xnames, drop = FALSE])
  
  # linear predictor for the NEW data
  lp <- as.vector(predict(object, newdata = X, type = "lp"))
  
  # grab the baseline cumulative hazard that we stored in the model;
  # build a step function so we can evaluate it at arbitrary 'times'
  H0_step <- stats::stepfun(object$bh_time,
                            c(0, object$bh_hazard))
  
  H0_vec  <- H0_step(times)          # baseline cumul. haz. at all eval times
  
  # matrix:  nrow = n_obs,  ncol = length(times)
  surv_mat <- exp(-outer(exp(lp), H0_vec, "*"))
  
  risk_mat <- 1 - surv_mat           # Score() expects RISKS
  return(risk_mat)
}

feature_data <- reactive({
  # Only compute this once when needed
  req(rsf_model)
  
  # Create a named vector for quick lookup
  gene_names <- character(length(lookup$safe))
  names(gene_names) <- lookup$safe
  
  # Process all gene names at once (much faster than individual lookups)
  for(i in seq_along(lookup$safe)) {
    original_id <- lookup$raw[i]
    if(original_id %in% rownames(GSE28735$featureData)) {
      gene_info <- GSE28735$featureData$gene_assignment[match(original_id, rownames(GSE28735$featureData))]
      gene_symbol <- extract_gene_symbol(gene_info)
      if(!is.na(gene_symbol)) {
        gene_names[lookup$safe[i]] <- paste0(gene_symbol, " (", lookup$safe[i], ")")
      } else {
        gene_names[lookup$safe[i]] <- lookup$safe[i]
      }
    } else {
      gene_names[lookup$safe[i]] <- lookup$safe[i]
    }
  }
  
  # Create the importance dataframe once
  importance_df <- data.frame(
    variable = names(rsf_model$importance),
    importance = rsf_model$importance,
    selected = names(rsf_model$importance) %in% selected_features,
    stringsAsFactors = FALSE
  )
  
  # Add display names using vectorized lookup (much faster than sapply)
  importance_df$display_name <- gene_names[importance_df$variable]
  
  return(importance_df)
})

withDefaultSpinner <- function(ui_element) {
  shinycssloaders::withSpinner(
    ui_element,
    type = 8,
    color = "#2E9FDF",
    size = 1
  )
}