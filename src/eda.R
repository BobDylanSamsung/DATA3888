# Applying summary function across samples
summary_values_eda <- summary(as.vector(gse$eMat))
summary_per_sample_eda <- as.data.frame(t(apply(gse$eMat, 2, summary)))
summary_per_sample_eda$SampleID <- rownames(summary_per_sample_eda)
summary_per_sample_eda <- summary_per_sample_eda %>%
  relocate(SampleID) %>%
  rename(
    Min = "Min.",
    Q1 = "1st Qu.",
    Median = "Median",
    Mean = "Mean",
    Q3 = "3rd Qu.",
    Max = "Max."
  )

summary_long_eda <- summary_per_sample_eda %>%
  pivot_longer(cols = c("Min", "Q1", "Median", "Mean", "Q3", "Max"), names_to = "Statistic", values_to = "Value")

# PCA Setup
gse$phenoData$Tissue <- factor(ifelse(grepl("non-tumor|normal", gse$phenoData$`source_name_ch1`, ignore.case = TRUE), "Normal", "Tumor"))
pca_all <- prcomp(t(gse$eMat), scale. = TRUE)

# Heatmap Data Preparation
var_genes_eda <- apply(gse$eMat, 1, var)
top100_genes_eda <- names(sort(var_genes_eda, decreasing = TRUE)[1:100])
annotation_col_eda <- data.frame(Tissue = gse$phenoData$Tissue)
rownames(annotation_col_eda) <- rownames(gse$phenoData)