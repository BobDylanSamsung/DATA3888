# Prepare training and test datasets
X <- t(gse$eMat) %>% as.matrix()
y <- as.numeric(gse$phenoData$is_dead)

train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_idx, ]
y_train <- y[train_idx]
X_test  <- X[-train_idx, ]
y_test  <- y[-train_idx]

# Fit model using cross-validation
cvfit_binom <- cv.glmnet(x = X_train, y = y_train, family = "binomial", alpha = 1, nfolds = 5)

# Save lambda value
best_lambda <- cvfit_binom$lambda.min

# Train final model with the optimal lambda
final_model_binom <- glmnet(x = X_train, y = y_train, family = "binomial", alpha = 1, lambda = best_lambda)
coef_binom <- coef(final_model_binom)
nz_idx <- which(coef_binom != 0)
df_genes <- data.frame(
  gene_id = rownames(coef_binom)[nz_idx],
  coef    = as.numeric(coef_binom[nz_idx])
)

# Predict on the test set
pred_test <- predict(final_model_binom, newx = X_test, type = "class")
pred_factor <- factor(pred_test, levels = c(0,1), labels = c("Alive","Dead"))
actual_factor <- factor(y_test, levels = c(0,1), labels = c("Alive","Dead"))

cm_results <- confusionMatrix(pred_factor, actual_factor, positive = "Dead")