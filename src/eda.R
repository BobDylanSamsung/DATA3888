# Combine pheno data
eda_pheno <- bind_rows(GSE28735$phenoData, GSE62452$phenoData)
eda_matrix <- combined_matrix[, colnames(combined_matrix) %in% rownames(eda_pheno)]
stopifnot(all(colnames(eda_matrix) == rownames(eda_pheno)))

gse_eda <- list(
  featureData = GSE28735$featureData,
  phenoData = eda_pheno,
  eMat = eda_matrix
)
summary_values_eda <- summary(as.vector(gse_eda$eMat))
# hist(gse_eda$eMat, breaks = 100, main = "All Samples Expression Histogram", xlab = "Expression", col = "lightblue")
summary_per_sample_eda <- as.data.frame(t(apply(gse_eda$eMat, 2, summary)))
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

kable(head(summary_per_sample_eda), caption = "Summary Statistics (All Samples)")
summary_long_eda <- summary_per_sample_eda %>%
  pivot_longer(cols = -SampleID, names_to = "Statistic", values_to = "Value")

#ggplot(summary_long_eda, aes(x = Statistic, y = Value)) +
  #geom_boxplot(fill = "skyblue") +
  #facet_wrap(~ Statistic, scales = "free") +
  #ggtitle("Summary Statistics Boxplot (All Samples)") +
  #theme_minimal()

#boxplot(gse_eda$eMat, outline = FALSE, las = 2, main = "Expression per Sample (All)", col = "gray90")
gse_eda$phenoData$Tissue <- ifelse(grepl("non-tumor|normal", gse_eda$phenoData$`source_name_ch1`, ignore.case = TRUE),
                                   "Normal", "Tumor")
gse_eda$phenoData$Tissue <- factor(gse_eda$phenoData$Tissue)

pca_all <- prcomp(t(gse_eda$eMat), scale. = TRUE)

autoplot(pca_all, data = gse_eda$phenoData, colour = "Tissue") +
  ggtitle("PCA: Tumor vs Normal Samples") +
  theme_minimal()
var_genes_eda <- apply(gse_eda$eMat, 1, var)
top100_genes_eda <- names(sort(var_genes_eda, decreasing = TRUE)[1:100])

annotation_col_eda <- data.frame(Tissue = gse_eda$phenoData$Tissue)
rownames(annotation_col_eda) <- rownames(gse_eda$phenoData)

pheatmap(gse_eda$eMat[top100_genes_eda, ],
         scale = "row",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = annotation_col_eda,
         main = "Top 100 Variable Genes - All Samples")