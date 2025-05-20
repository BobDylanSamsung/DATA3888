# Load required libraries
library(caret)
library(CoxBoost)
library(dplyr)
library(DT)
library(GEOquery)
library(ggplot2)
library(ggfortify)
library(glmnet)
library(Hmisc)
library(knitr)
library(limma)
library(pheatmap)
library(randomForestSRC)
library(pec)
library(pROC)
library(prodlim)
library(RColorBrewer)
library(riskRegression)
library(rms)
library(reshape2)
library(shiny)
library(shinydashboard)
library(survival)
library(survminer)
library(tidyverse)
library(timeROC)

set.seed(3888)

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
