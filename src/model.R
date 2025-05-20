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

message("[INFO]: Selecting Features")

# Run Random Survival Forest for feature selection
rsf_model <- rfsrc(Surv(time, is_dead) ~ ., data = train_df, 
                   ntree = 1000, 
                   importance = TRUE)

# Get variable importance by minimal depth
var_imp <- var.select(rsf_model, method = "md", verbose = FALSE)

top_by_md <- var_imp$md.obj$topvars          # may be character(0)

if (length(top_by_md) == 0) {                # fall-back: use VIMP
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

# STANDARDIZATION STEP
X_train_selected_std <- scale(X_train_selected)
X_test_selected_std <- scale(X_test_selected, 
                             center = attr(X_train_selected_std, "scaled:center"),
                             scale = attr(X_train_selected_std, "scaled:scale"))

# Set up CoxBoost - note it requires time and status separately
train_time <- train_df$time
train_is_dead <- train_df$is_dead  
test_time <- test_df$time
test_is_dead <- test_df$is_dead    

# choose a reasonable penalty
penalty_val <- 4

message("[INFO]: Calculating optimal boosting steps")
# Cross validate to find optimal number of boosting steps
cv_res <- cv.CoxBoost(
  time      = train_time,
  status    = train_is_dead,
  x         = X_train_selected_std,
  maxstepno = 100,
  K         = 10,
  penalty   = penalty_val,
  unpen.index = NULL
)

# extract optimal steps
optimal_steps <- cv_res$optimal.step
message("[INFO]: Optimal boosting steps = ", optimal_steps)
message("[INFO]: Fitting model")
# Fit final CoxBoost model with the same penalty
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
          at.step = optimal_steps)   # <- returns a vector
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

# Generate confusion matrix
cm_caret <- confusionMatrix(
  data = risk_groups,
  reference = actual_is_dead,
  positive = "High"
)