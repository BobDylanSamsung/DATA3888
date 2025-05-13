# Function to extract gene symbol from gene assignment
extract_gene_symbol <- function(gene_assignment) {
  # Split the string by //
  parts <- strsplit(gene_assignment, " // ")[[1]]
  
  # Return the second element (gene symbol)
  if (length(parts) >= 2) {
    return(trimws(parts[2]))
  } else {
    return(NA)
  }
}

# Function to extract data from GEO series
extract_gse_data <- function(gse) {
  return(list(
    featureData = fData(gse),
    phenoData = pData(gse),
    eMat = exprs(gse)
  ))
}

# Robust function to download GEO data with error handling
get_geo_data <- function(gse_id, destdir = tempdir(), retries = 5) {
  message(paste("[INFO]: Attempting to load", gse_id))
  
  for (i in 1:retries) {
    tryCatch({
      # Try to download the data
      gse_list <- getGEO(gse_id, destdir = destdir, GSEMatrix = TRUE, AnnotGPL = FALSE)
      
      # Extract the first element (series matrix)
      gse <- gse_list[[1]]
      
      message(paste("[INFO]: Successfully downloaded", gse_id))
      return(gse)
      
    }, error = function(e) {
      message(paste("[INFO]: Attempt", i, "failed for", gse_id, ":", e$message))
      
      if (i < retries) {
        message("Retrying in 3 seconds...")
        Sys.sleep(3)
      } else {
        stop(paste("Failed to download", gse_id, "after", retries, "attempts"))
      }
    })
  }
}

download_dir <- file.path(getwd(), "GEO_data")
if (!dir.exists(download_dir)) {
  dir.create(download_dir)
}

# Load datasets
tryCatch({
  GSE28735 <- get_geo_data("GSE28735", destdir = download_dir) %>% extract_gse_data()
  GSE62452 <- get_geo_data("GSE62452", destdir = download_dir) %>% extract_gse_data()
}, error = function(e) {
  message(paste("Error loading additional datasets:", e$message))
})

common_genes <- intersect(rownames(GSE28735$eMat), rownames(GSE62452$eMat))
expr_combined <- cbind(
  GSE28735$eMat[common_genes, ],
  GSE62452$eMat[common_genes, ]
)

# Phenotype processing
pheno287 <- GSE28735$phenoData %>%
  mutate(is_dead = as.numeric(replace(`cancer_death:ch1`, `cancer_death:ch1` == "na", NA))) %>%
  filter(!grepl("non-tumor tissue", source_name_ch1) & !is.na(is_dead)) %>%
  rename(time = `survival_month:ch1`)

pheno624 <- GSE62452$phenoData %>%
  mutate(is_dead = as.numeric(replace(`survival status:ch1`, `survival status:ch1` %in% c("na","?"), NA))) %>%
  filter(grepl("Pancreatic tumor", `tissue:ch1`) & !is.na(is_dead)) %>%
  rename(time = `survival months:ch1`)

pheno <- bind_rows(pheno287, pheno624)
expr_mat <- expr_combined[, colnames(expr_combined) %in% rownames(pheno)]

# Remove a row for each is_dead for later use
dead_0_sample <- rownames(pheno)[pheno$is_dead == 0][4]
dead_1_sample <- rownames(pheno)[pheno$is_dead == 1][4]
expr_0 <- expr_mat[, dead_0_sample, drop = FALSE]
expr_1 <- expr_mat[, dead_1_sample, drop = FALSE]
write.csv(expr_0, file = "predict_test_dead0.csv")
write.csv(expr_1, file = "predict_test_dead1.csv")

exclude_samples <- c(dead_0_sample, dead_1_sample)
pheno <- pheno[!rownames(pheno) %in% exclude_samples, ]
expr_mat <- expr_mat[, !colnames(expr_mat) %in% exclude_samples]

subset_eMat_GSE28735 <- GSE28735$eMat[common_genes, ]
subset_eMat_GSE62452 <- GSE62452$eMat[common_genes, ]
combined_matrix <- cbind(subset_eMat_GSE28735, subset_eMat_GSE62452)
combined_matrix <- combined_matrix[, colnames(combined_matrix) %in% rownames(pheno)]
stopifnot(all(colnames(combined_matrix) == rownames(pheno)))

gse <- list(
  featureData = GSE28735$featureData,
  phenoData = pheno,
  eMat = combined_matrix
)
