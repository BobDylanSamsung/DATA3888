# Construct survival object with survival time and event status
surv_obj <- with(gse$phenoData, Surv(as.numeric(months_survived), as.numeric(is_dead)))

# Prepare feature matrix
X <- t(gse$eMat) %>% as.matrix()

# Split into training and testing sets
train_idx <- createDataPartition(gse$phenoData$is_dead, p = 0.8, list = FALSE)
X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]

# Survival objects for training/testing sets
train_surv <- surv_obj[train_idx]
test_surv <- surv_obj[-train_idx]

# Fit penalized Cox model using cross-validation
cvfit_cox <- cv.glmnet(x = X_train, y = train_surv, family = "cox", alpha = 1, nfolds = 5)

# Select the optimal lambda
best_lambda_cox <- cvfit_cox$lambda.min

# Fit final Cox model with optimal lambda
final_model_cox <- glmnet(x = X_train, y = train_surv, family = "cox", alpha = 1, lambda = best_lambda_cox)

# Extract coefficients from the final model
coef_cox <- coef(final_model_cox)
nz_idx_cox <- which(coef_cox != 0)
df_genes_cox <- data.frame(
  gene_id = colnames(X_train)[nz_idx_cox],
  coef    = as.numeric(coef_cox[nz_idx_cox])
)

# Use the model to predict risk scores on the test set (linear predictors)
risk_scores_test <- predict(final_model_cox, newx = X_test, type = "link")

# Split phenotype data into training and testing sets
train_pheno <- gse$phenoData[train_idx, ]
test_pheno <- gse$phenoData[-train_idx, ]

# Compute risk scores for each test sample, using precomputed risk scores if applicable
test_pheno$predicted_risk <- ifelse(risk_scores_test > median(risk_scores_test), "High", "Low")
test_pheno$predicted_risk <- factor(test_pheno$predicted_risk, levels = c("Low", "High"))

# Construct actual labels for evaluation (risk groups using median)
test_pheno$actual_status <- factor(
  ifelse(
    test_pheno$is_dead == 1,
    "High",
    "Low"
  ),
  levels = c("Low", "High") 
)

# Calculate C-index for model performance evaluation on test set
c_index <- rcorr.cens(-risk_scores_test, test_surv)[1]

# Generate confusion matrix using caret
cm_caret <- confusionMatrix(
  data = test_pheno$predicted_risk,
  reference = test_pheno$actual_status,
  positive = "High"
)