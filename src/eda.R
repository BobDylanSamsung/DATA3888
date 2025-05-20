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

# Heatmap Data Preparation - FIX: Increase to 75 genes
var_genes_eda <- apply(gse$eMat, 1, var)
top75_genes_eda <- names(sort(var_genes_eda, decreasing = TRUE)[1:75])
annotation_col_eda <- data.frame(Survival = gse$phenoData$Survival_Status)
rownames(annotation_col_eda) <- rownames(gse$phenoData)