message("[INFO]: Constructing survival object...")
surv_obj <- with(pheno, Surv(as.numeric(time), is_dead))

message("[INFO]: Preparing expression matrix...")
X <- t(expr_mat) %>% as.matrix()

# Lasso-Cox for feature selection
message("[INFO]: Performing cross-validation for Lasso-Cox model...")
set.seed(3888)
cvfit <- cv.glmnet(X, surv_obj, family = "cox", alpha = 1)

message("[INFO]: Extracting best lambda from cross-validation...")
best_lambda <- cvfit$lambda.min

message("[INFO]: Building Lasso-Cox model with best lambda...")
lasso_model <- glmnet(X, surv_obj, family = "cox", alpha = 1, lambda = best_lambda)

message("[INFO]: Extracting non-zero coefficients for selected genes...")
gene_coefs <- data.frame(
  gene = colnames(X)[coef(lasso_model)@i + 1],
  coef = coef(lasso_model)@x
)
selected_genes <- gene_coefs$gene

message("[INFO]: Preparing data frame for CoxPH model refitting...")
train_df <- data.frame(
  time   = as.numeric(pheno$time),
  status = as.numeric(pheno$is_dead),
  as.data.frame(X[, selected_genes, drop = FALSE])
)
colnames(train_df)[3:ncol(train_df)] <- selected_genes

message("[INFO]: Refitting the CoxPH model...")
coxph_model <- coxph(Surv(time, status) ~ ., data = train_df)

message("[INFO]: Calculating risk scores and determining cut-off...")
train_lp <- predict(coxph_model, newdata = train_df, type = "lp")
cutoff <- median(train_lp, na.rm = TRUE)