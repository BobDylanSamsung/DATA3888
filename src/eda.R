# Applying summary function across samples
summary_values_eda <- summary(as.vector(gse$eMat))
summary_per_sample_eda <- as.data.frame(t(apply(gse$eMat, 2, summary)))

# Set SampleID as a column (not using rownames to avoid duplication)
summary_per_sample_eda$SampleID <- rownames(summary_per_sample_eda)
rownames(summary_per_sample_eda) <- NULL  # Remove rownames to prevent duplication

summary_per_sample_eda <- summary_per_sample_eda %>%
  rename(
    Min = "Min.",
    Q1 = "1st Qu.",
    Median = "Median",
    Mean = "Mean",
    Q3 = "3rd Qu.",
    Max = "Max."
  ) %>%
  
  select(SampleID, Min, Q1, Median, Mean, Q3, Max)

summary_long_eda <- summary_per_sample_eda %>%
  pivot_longer(cols = c("Min", "Q1", "Median", "Mean", "Q3", "Max"), names_to = "Statistic", values_to = "Value")

# PCA Setup - using survival status instead of tissue
gse$phenoData$Survival_Status <- factor(
  ifelse(gse$phenoData$is_dead == 0, "Alive", "Deceased"),
  levels = c("Alive", "Deceased")
)

pca_all <- prcomp(t(gse$eMat), scale. = TRUE)

# Heatmap Data Preparation 
var_genes_eda <- apply(gse$eMat, 1, var)
top75_genes_eda <- names(sort(var_genes_eda, decreasing = TRUE)[1:75])
annotation_col_eda <- data.frame(Survival = gse$phenoData$Survival_Status)
rownames(annotation_col_eda) <- rownames(gse$phenoData)

# Prepare data for survival time distribution
prepare_survival_time_data <- function() {
  # Create a dataframe with survival time and status
  survival_data <- data.frame(
    SampleID = rownames(gse$phenoData),
    SurvivalTime = gse$phenoData$months_survived, 
    Status = factor(
      ifelse(gse$phenoData$is_dead == 0, "Alive", "Deceased"),
      levels = c("Alive", "Deceased")
    )
  )
  return(survival_data)
}

# Prepare data for volcano plot - differential expression between alive/deceased
prepare_volcano_data <- function() {
  # Create design matrix
  design <- model.matrix(~gse$phenoData$is_dead)
  colnames(design) <- c("Intercept", "Deceased_vs_Alive")
  
  # Fit linear model
  fit <- limma::lmFit(gse$eMat, design)
  fit <- limma::eBayes(fit)
  
  # Get results
  results <- limma::topTable(fit, coef = "Deceased_vs_Alive", number = Inf)
  
  # Add gene names
  results$GeneID <- rownames(results)
  
  # Flag significant genes (adjust threshold as needed)
  results$Significant <- ifelse(results$adj.P.Val < 0.05, "Yes", "No")
  results$Label <- ifelse(results$Significant == "Yes" & abs(results$logFC) > 1, results$GeneID, "")
  
  return(results)
}

survival_time_data <- prepare_survival_time_data()
volcano_data <- prepare_volcano_data()