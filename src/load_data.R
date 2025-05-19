# --- Helper Functions ---

# Extract gene symbol from gene_assignment string
extract_gene_symbol <- function(gene_assignment) {
  parts <- strsplit(gene_assignment, " // ")[[1]]
  if (length(parts) >= 2) return(trimws(parts[2])) else return(NA)
}

# Extract useful data from GEO object
extract_gse_data <- function(gse) {
  list(
    featureData = fData(gse),
    phenoData = pData(gse),
    eMat = exprs(gse)
  )
}

# Robust GEO downloader
get_geo_data <- function(gse_id, destdir = tempdir(), retries = 5) {
  message(paste("[INFO]: Attempting to load", gse_id))
  for (i in 1:retries) {
    tryCatch({
      gse_list <- getGEO(gse_id, destdir = destdir, GSEMatrix = TRUE, AnnotGPL = FALSE)
      gse <- gse_list[[1]]
      message(paste("[INFO]: Successfully downloaded", gse_id))
      return(gse)
    }, error = function(e) {
      message(paste("[INFO]: Attempt", i, "failed for", gse_id, ":", e$message))
      if (i < retries) Sys.sleep(3) else stop(paste("Failed to download", gse_id))
    })
  }
}

# --- Step 1: Load GEO datasets ---
download_dir <- file.path(getwd(), "GEO_data")
if (!dir.exists(download_dir)) dir.create(download_dir)

GSE28735 <- get_geo_data("GSE28735", destdir = download_dir) %>% extract_gse_data()
GSE62452 <- get_geo_data("GSE62452", destdir = download_dir) %>% extract_gse_data()

# --- Step 2: Clean phenotype data ---

pheno287 <- GSE28735$phenoData %>%
  mutate(
    is_dead = as.numeric(replace(`cancer_death:ch1`, `cancer_death:ch1` == "na", NA)),
    months_survived = as.numeric(replace(`survival_month:ch1`, `survival_month:ch1` == "na", NA))
  ) %>%
  filter(!grepl("non-tumor tissue", source_name_ch1) & !is.na(is_dead))

pheno624 <- GSE62452$phenoData %>%
  mutate(
    is_dead = as.numeric(replace(`survival status:ch1`, `survival status:ch1` %in% c("na","?"), NA)),
    months_survived = as.numeric(replace(`survival months:ch1`, `survival months:ch1` %in% c("na","?"), NA))
  ) %>%
  filter(grepl("Pancreatic tumor", `tissue:ch1`) & !is.na(is_dead))

# Combine phenotype
pheno <- bind_rows(pheno287, pheno624)

# --- Step 3: Combine expression data ---

common_genes <- intersect(rownames(GSE28735$eMat), rownames(GSE62452$eMat))
subset_eMat_GSE28735 <- GSE28735$eMat[common_genes, ]
subset_eMat_GSE62452 <- GSE62452$eMat[common_genes, ]
combined_matrix <- cbind(subset_eMat_GSE28735, subset_eMat_GSE62452)

# --- Step 4: Save validation samples and remove them from training data ---

selected_ids <- c("GSM1527137", "GSM1527230", "GSM711944", "GSM711984")

# Save selected samples to CSV
outdir <- file.path(getwd(), "patient_data")
if (!dir.exists(outdir)) dir.create(outdir)

for (sample_id in selected_ids) {
  expr_df <- data.frame(
    gene_id    = rownames(combined_matrix),
    expression = combined_matrix[, sample_id]
  )
  ph <- pheno[sample_id, ]
  expr_df$is_dead         <- ph$is_dead
  expr_df$months_survived <- ph$months_survived
  
  write.csv(expr_df, file = file.path(outdir, paste0(sample_id,"is_dead_", ph$is_dead, ".csv")),
            row.names = FALSE, quote = FALSE)
}

# Remove validation samples
pheno <- pheno[!rownames(pheno) %in% selected_ids, ]
combined_matrix <- combined_matrix[, colnames(combined_matrix) %in% rownames(pheno)]

# Final check
stopifnot(all(colnames(combined_matrix) == rownames(pheno)))

# Final object
gse <- list(
  featureData = GSE28735$featureData,
  phenoData = pheno,
  eMat = combined_matrix
)
