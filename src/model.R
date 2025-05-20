# Construct survival object with survival time and event status
surv_obj <- with(gse$phenoData, Surv(as.numeric(months_survived), as.numeric(is_dead)))

# Prepare feature matrix
X <- t(gse$eMat) %>% as.matrix()
train_idx <- createDataPartition(factor(gse$phenoData$is_dead), p = 0.80, list = FALSE)

safe_names <- make.names(colnames(X), unique = TRUE)
colnames(X) <- safe_names
lookup <- data.frame(
  safe = safe_names,              # length 28869
  raw  = rownames(gse$eMat),      # length 28869
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
  as.data.frame(X_test , check.names = FALSE)
)

message("[INFO]: Starting feature selection with Random Survival Forest...")

# Run Random Survival Forest for feature selection
message("[INFO]: Growing 1000 trees for RSF model...")
rsf_model <- rfsrc(Surv(time, is_dead) ~ ., data = train_df, 
                   ntree = 1000, 
                   importance = TRUE)
message("[INFO]: RSF model completed. Calculating variable importance...")

# Get variable importance by minimal depth
var_imp <- var.select(rsf_model, method = "md", verbose = FALSE)

top_by_md <- var_imp$md.obj$topvars          # may be character(0)

if (length(top_by_md) == 0) {                # fall-back: use VIMP
  message("[INFO]: No variables selected by minimal depth. Using VIMP instead...")
  vimp_sorted <- sort(rsf_model$importance,
                      decreasing = TRUE,
                      na.last   = TRUE)
  top_by_md <- names(vimp_sorted)
}

n_features_to_select <- 15
selected_features <- top_by_md[1:min(n_features_to_select,
                                     length(top_by_md))]

if (length(selected_features) == 0)
  stop("[ERROR]: variable-selection produced 0 features")

message("[INFO]: Number of features selected: ", length(selected_features))

# Prepare data for CoxBoost with selected features only
X_train_selected <- X_train[, selected_features, drop = FALSE]
X_test_selected <- X_test[, selected_features, drop = FALSE]

# Standardize the data
X_train_selected_std <- scale(X_train_selected)
X_test_selected_std <- scale(X_test_selected, 
                             center = attr(X_train_selected_std, "scaled:center"),
                             scale = attr(X_train_selected_std, "scaled:scale"))

# Set up CoxBoost - note it requires time and status separately
train_time <- train_df$time
train_is_dead <- train_df$is_dead  
test_time <- test_df$time
test_is_dead <- test_df$is_dead    

# Hyperparameter tuning for CoxBoost
message("[INFO]: Starting hyperparameter tuning for CoxBoost...")

# Define a grid of penalty values to try
penalty_grid <- seq(from = 0, to = 50, by = 5)
max_steps <- 100

# Store results
tuning_results <- data.frame(penalty = numeric(), 
                             optimal_steps = numeric(), 
                             c_index = numeric(),
                             loglik = numeric(),
                             stringsAsFactors = FALSE)

# Perform grid search
for (penalty_val in penalty_grid) {
  message("[INFO]: Testing penalty value: ", penalty_val)
  
  # Cross validate to find optimal number of boosting steps for this penalty
  cv_res <- cv.CoxBoost(
    time      = train_time,
    status    = train_is_dead,
    x         = X_train_selected_std,
    maxstepno = max_steps,
    K         = 5,
    penalty   = penalty_val,
    unpen.index = NULL
  )
  
  # extract optimal steps
  optimal_steps <- cv_res$optimal.step
  
  message("[INFO]: Penalty ", penalty_val, " - Optimal boosting steps = ", optimal_steps)
  
  # Get CV likelihood based on the CV result (with check for NULL or missing)
  cv_loglik <- if(is.null(cv_res$cvlp) || length(cv_res$cvlp) == 0) NA else cv_res$cvlp
  
  # Calculate C-index for this model configuration
  temp_model <- CoxBoost(
    time      = train_time,
    status    = train_is_dead,
    x         = X_train_selected_std,
    stepno    = optimal_steps,
    penalty   = penalty_val,
    unpen.index = NULL
  )
  
  # Create validation set
  val_idx <- sample(length(train_time), round(length(train_time) * 0.2))
  val_preds <- predict(temp_model, 
                       newdata = X_train_selected_std[val_idx,], 
                       type = "lp", 
                       at.step = optimal_steps)
  
  # Calculate C-index with error handling
  c_index_val <- tryCatch({
    rcorr.cens(-val_preds, Surv(train_time[val_idx], train_is_dead[val_idx]))[1]
  }, error = function(e) {
    message("[WARNING]: Could not calculate C-index: ", e$message)
    return(NA)
  })
  
  # Store results with error handling for all variables
  tuning_results <- rbind(tuning_results, 
                          data.frame(penalty = penalty_val, 
                                     optimal_steps = optimal_steps, 
                                     c_index = c_index_val,
                                     loglik = cv_loglik,
                                     stringsAsFactors = FALSE))
  
  if(!is.na(c_index_val)) {
    message("[INFO]: Penalty ", penalty_val, " validation C-index: ", round(c_index_val, 4))
  } else {
    message("[INFO]: Penalty ", penalty_val, " - Could not calculate validation C-index")
  }
}

message("[INFO]: Hyperparameter tuning completed.")
print(tuning_results)

# Find best hyperparameters
best_params <- tuning_results[which.max(tuning_results$c_index), ]
message("[INFO]: Best hyperparameters - Penalty: ", best_params$penalty, 
        ", Steps: ", best_params$optimal_steps, 
        ", C-index: ", round(best_params$c_index, 4))

# Use best hyperparameters for final model
penalty_val <- best_params$penalty
optimal_steps <- best_params$optimal_steps

message("[INFO]: Fitting final model with best hyperparameters")
# Fit final CoxBoost model with best hyperparameters
coxboost_model <- CoxBoost(
  time      = train_time,
  status    = train_is_dead,
  x         = X_train_selected_std,
  stepno    = optimal_steps,
  penalty   = penalty_val,
  unpen.index = NULL
)

# linear predictor for the training data
lp_train <- as.vector(
  predict(coxboost_model,
          newdata = X_train_selected_std,
          type    = "lp",
          at.step = optimal_steps)   # returns a vector
)

# data frame that will be fed to coxph()
df_train <- data.frame(
  time   = train_time,
  status = train_is_dead,
  lp     = lp_train       
)

# "null" Cox model 
coxph_base <- coxph(
  Surv(time, status) ~ offset(lp),   
  data   = df_train,
  method = "breslow"
)

bh <- basehaz(coxph_base, centered = FALSE)

coxboost_model$bh_time   <- bh$time
coxboost_model$bh_hazard <- bh$hazard
coxboost_model$stepno <- optimal_steps

# Predict on test set
risk_scores_test <- as.vector(
  predict(coxboost_model,
          newdata = X_test_selected_std,
          type    = "lp",
          at.step = optimal_steps)
)

# Split test data into high/low risk groups using median cutoff
risk_groups <- ifelse(risk_scores_test > median(risk_scores_test), "High", "Low")
risk_groups <- factor(risk_groups, levels = c("Low", "High"))

# Create actual status for comparison
actual_is_dead <- factor(
  ifelse(test_is_dead == 1, "High", "Low"),
  levels = c("Low", "High")
)

# Calculate C-index
c_index <- rcorr.cens(-risk_scores_test, Surv(test_time, test_is_dead))[1]
message("[INFO]: Final model C-index on test data: ", round(c_index, 4))

# Generate confusion matrix
cm_caret <- confusionMatrix(
  data = risk_groups,
  reference = actual_is_dead,
  positive = "High"
)

