# Script to find the optimal CoxBoost parameters for survival analysis
source("src/global.R")
source("src/load_data.R")

# Construct survival object
surv_obj <- with(gse$phenoData, Surv(as.numeric(months_survived), as.numeric(is_dead)))

# Prepare feature matrix
X <- t(gse$eMat) %>% as.matrix()
train_idx <- createDataPartition(factor(gse$phenoData$is_dead), p = 0.80, list = FALSE)

safe_names <- make.names(colnames(X), unique = TRUE)
colnames(X) <- safe_names
lookup <- data.frame(
  safe = safe_names,
  raw  = rownames(gse$eMat),
  stringsAsFactors = FALSE
)

# Split into training and testing sets
X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]

# Survival objects for training/testing sets
train_surv <- surv_obj[train_idx]
test_surv <- surv_obj[-train_idx]

# Create dataframes for training and testing
train_df <- data.frame(
  time    = as.numeric(gse$phenoData$months_survived[train_idx]),
  is_dead = as.numeric(gse$phenoData$is_dead[train_idx]),
  as.data.frame(X_train, check.names = FALSE)
)

test_df  <- data.frame(
  time    = as.numeric(gse$phenoData$months_survived[-train_idx]),
  is_dead = as.numeric(gse$phenoData$is_dead[-train_idx]),
  as.data.frame(X_test, check.names = FALSE)
)

# Feature selection with Random Survival Forest
message("[INFO]: Growing 1000 trees for RSF model...")
rsf_model <- rfsrc(Surv(time, is_dead) ~ ., data = train_df, 
                   ntree = 1000, 
                   importance = TRUE)

# Get variable importance by minimal depth
var_imp <- var.select(rsf_model, method = "md", verbose = FALSE)
top_by_md <- var_imp$md.obj$topvars

if (length(top_by_md) == 0) {
  message("[INFO]: No variables selected by minimal depth. Using VIMP instead...")
  vimp_sorted <- sort(rsf_model$importance, decreasing = TRUE, na.last = TRUE)
  top_by_md <- names(vimp_sorted)
}

n_features_to_select <- 15
selected_features <- top_by_md[1:min(n_features_to_select, length(top_by_md))]

# Prepare data for CoxBoost with selected features only
X_train_selected <- X_train[, selected_features, drop = FALSE]
X_test_selected <- X_test[, selected_features, drop = FALSE]

# Standardize the data
X_train_selected_std <- scale(X_train_selected)
X_test_selected_std <- scale(X_test_selected, 
                             center = attr(X_train_selected_std, "scaled:center"),
                             scale = attr(X_train_selected_std, "scaled:scale"))

# Set up CoxBoost variables
train_time <- train_df$time
train_is_dead <- train_df$is_dead  
test_time <- test_df$time
test_is_dead <- test_df$is_dead

# Define a helper function to evaluate CoxBoost model with Brier score
evaluate_coxboost <- function(penalty_val, optimal_steps, X_train_std, X_test_std, 
                              train_time, train_is_dead, test_time, test_is_dead, selected_features) {
  # Fit CoxBoost model
  coxboost_model <- CoxBoost(
    time = train_time,
    status = train_is_dead,
    x = X_train_std,
    stepno = optimal_steps,
    penalty = penalty_val,
    unpen.index = NULL
  )
  
  # Add baseline hazard to model for Brier score calculation
  lp_train <- as.vector(
    predict(coxboost_model,
            newdata = X_train_std,
            type = "lp",
            at.step = optimal_steps)
  )
  
  df_train <- data.frame(
    time = train_time,
    status = train_is_dead,
    lp = lp_train       
  )
  
  coxph_base <- coxph(
    Surv(time, status) ~ offset(lp),   
    data = df_train,
    method = "breslow"
  )
  
  bh <- basehaz(coxph_base, centered = FALSE)
  
  coxboost_model$bh_time <- bh$time
  coxboost_model$bh_hazard <- bh$hazard
  coxboost_model$stepno <- optimal_steps
  coxboost_model$xnames <- selected_features
  
  # Calculate C-index
  risk_scores_test <- as.vector(
    predict(coxboost_model,
            newdata = X_test_std,
            type = "lp",
            at.step = optimal_steps)
  )
  
  c_index <- rcorr.cens(-risk_scores_test, Surv(test_time, test_is_dead))[1]
  
  # Calculate Brier score
  eval_times <- seq(6, floor(max(test_time)), by = 6)
  eval_times <- eval_times[sapply(eval_times, function(t)
    any(test_is_dead == 1 & test_time <= t) && any(test_time >= t))]
  
  test_score_df <- data.frame(
    time = test_time,
    is_dead = test_is_dead
  )
  
  # Add feature columns
  for (feat in selected_features) {
    test_score_df[[feat]] <- X_test_std[, feat]
  }
  
  # Calculate Brier score with error handling
  score_res <- tryCatch({
    riskRegression::Score(
      object = list(CoxBoost = coxboost_model),
      formula = Surv(time, is_dead) ~ 1,
      data = test_score_df,
      times = eval_times,
      metrics = "brier",
      cens.model = "km",
      split.method = "none",
      summary = "ibs"
    )
  }, error = function(e) {
    message("Error in Brier score calculation: ", e$message)
    return(NULL)
  })
  
  # Process results
  if (is.null(score_res) || length(eval_times) < 2) {
    ibs_value <- NA_real_
  } else {
    ibs_value <- score_res$Brier$score[score_res$Brier$score$model == "CoxBoost" & 
                                         score_res$Brier$score$times == max(score_res$Brier$score$times), "IBS"]
  }
  
  return(list(
    c_index = c_index,
    ibs_value = ibs_value,
    model = coxboost_model
  ))
}

# Hyperparameter tuning grid
message("[INFO]: Starting hyperparameter tuning for CoxBoost...")
penalty_grid <- seq(from = 0, to = 50, by = 5)
max_steps <- 100

# Store results
tuning_results <- data.frame(
  penalty = numeric(), 
  optimal_steps = numeric(), 
  cv_loglik = numeric(),
  c_index = numeric(),
  ibs = numeric(),
  stringsAsFactors = FALSE
)

# Perform grid search
for (penalty_val in penalty_grid) {
  message("[INFO]: Testing penalty value: ", penalty_val)
  
  # Cross validate to find optimal number of boosting steps
  cv_res <- cv.CoxBoost(
    time = train_time,
    status = train_is_dead,
    x = X_train_selected_std,
    maxstepno = max_steps,
    K = 5,
    penalty = penalty_val,
    unpen.index = NULL
  )
  
  # Extract optimal steps
  optimal_steps <- cv_res$optimal.step
  
  # Get cross-validated log-likelihood at optimal step
  cv_loglik <- cv_res$mean.logplik[optimal_steps + 1]  # +1 because index starts at 0
  
  # Evaluate model with C-index and Brier score
  eval_result <- evaluate_coxboost(
    penalty_val = penalty_val,
    optimal_steps = optimal_steps,
    X_train_std = X_train_selected_std,
    X_test_std = X_test_selected_std,
    train_time = train_time,
    train_is_dead = train_is_dead,
    test_time = test_time,
    test_is_dead = test_is_dead,
    selected_features = selected_features
  )
  
  # Store results
  tuning_results <- rbind(tuning_results, 
                          data.frame(
                            penalty = penalty_val, 
                            optimal_steps = optimal_steps, 
                            cv_loglik = cv_loglik,
                            c_index = eval_result$c_index,
                            ibs = eval_result$ibs_value,
                            stringsAsFactors = FALSE
                          ))
  
  message("[INFO]: Penalty ", penalty_val, 
          " - Optimal steps = ", optimal_steps,
          " - CV log-likelihood = ", round(cv_loglik, 4),
          " - C-index = ", round(eval_result$c_index, 4),
          " - IBS = ", round(eval_result$ibs_value, 4))
}

# Find best hyperparameters using multiple metrics
message("\n[INFO]: Hyperparameter tuning completed.")
print(tuning_results)

# Find best model by each metric
best_by_loglik <- tuning_results[which.max(tuning_results$cv_loglik), ]
best_by_cindex <- tuning_results[which.max(tuning_results$c_index), ]
best_by_ibs <- tuning_results[which.min(tuning_results$ibs), ]  # Lower IBS is better

message("\n[INFO]: Best by log-likelihood: Penalty = ", best_by_loglik$penalty, 
        ", Steps = ", best_by_loglik$optimal_steps,
        ", CV log-likelihood = ", round(best_by_loglik$cv_loglik, 4),
        ", C-index = ", round(best_by_loglik$c_index, 4),
        ", IBS = ", round(best_by_loglik$ibs, 4))

message("[INFO]: Best by C-index: Penalty = ", best_by_cindex$penalty, 
        ", Steps = ", best_by_cindex$optimal_steps,
        ", CV log-likelihood = ", round(best_by_cindex$cv_loglik, 4),
        ", C-index = ", round(best_by_cindex$c_index, 4),
        ", IBS = ", round(best_by_cindex$ibs, 4))

message("[INFO]: Best by IBS: Penalty = ", best_by_ibs$penalty, 
        ", Steps = ", best_by_ibs$optimal_steps,
        ", CV log-likelihood = ", round(best_by_ibs$cv_loglik, 4),
        ", C-index = ", round(best_by_ibs$c_index, 4),
        ", IBS = ", round(best_by_ibs$ibs, 4))

# Write best parameters to a text file for easy reference
writeLines(
  c("# CoxBoost Best Parameters",
    paste("# Generated on:", Sys.time()),
    "",
    "# Best parameters by log-likelihood:",
    paste("BEST_PENALTY_LOGLIK =", best_by_loglik$penalty),
    paste("BEST_STEPS_LOGLIK =", best_by_loglik$optimal_steps),
    "",
    "# Best parameters by C-index:",
    paste("BEST_PENALTY_CINDEX =", best_by_cindex$penalty),
    paste("BEST_STEPS_CINDEX =", best_by_cindex$optimal_steps),
    "",
    "# Best parameters by Integrated Brier Score:",
    paste("BEST_PENALTY_IBS =", best_by_ibs$penalty),
    paste("BEST_STEPS_IBS =", best_by_ibs$optimal_steps),
    "",
    "# Selected features:",
    paste(selected_features, collapse = ", ")
  ),
  "coxboost_best_params.txt"
)

message("[INFO]: Results saved to coxboost_tuning_results.RData and coxboost_best_params.txt")