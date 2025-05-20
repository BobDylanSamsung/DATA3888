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
