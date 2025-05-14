
# Pancreatic Cancer Survival Prediction Tool

This Shiny web application provides **personalized survival predictions** for pancreatic cancer patients based on gene expression data. Built on high-throughput GEO datasets and a **LASSO-Cox proportional hazards model**, the tool delivers interpretable risk assessments for precision medicine applications.

![screenshot](./assets/app_screenshot.png)

---

## 🚀 Features

- 📊 **Exploratory Data Analysis (EDA)**  
  Visualizations and summaries of gene expression patterns across datasets.

- 🧠 **Model Training & Evaluation**  
  LASSO-Cox regression with cross-validation and survival performance metrics.

- 🧬 **Patient Risk Prediction**  
  Upload a CSV of gene expression for a new patient to receive:
  - Risk Score
  - Risk Group (High/Low)
  - Hazard Ratio
  - Estimated Median Survival
  - Top Contributing Genes
  - Risk Percentile Visualization

---

## 📁 Project Structure

```
.
├── data/                 # Processed gene expression and phenotype data
├── models/               # Trained Lasso-Cox model object
├── src/
│   ├── ui/               # Modular Shiny UI components
│   │   ├── ui.R
│   │   ├── home.R
│   ├── server.R          # Shiny server logic
│   └── global.R          # Global setup (model loading, constants, etc.)
├── assets/               # image assets
├── app.R                 # Main Shiny app entry point
└── README.md             # This file
```

---

## 📥 Input Format

To generate a prediction, upload a `.csv` file with the following structure:

```csv
Gene,Expression
ID1,2.43
ID2,4.55
...
```

- **Gene:** HGNC gene ID (must match those used in training)
- **Expression:** Normalized expression value

---

## ⚙️ Installation

1. Clone the repository:

```bash
git clone https://github.com/BobDylanSamsung/DATA3888.git
cd DATA3888
```

2. Install dependencies:

Run the following code to install all required packages:

```r
install.packages(c(
  "caret", "dplyr", "GEOquery", "ggplot2", "ggfortify", 
  "glmnet", "Hmisc", "knitr", "limma", "pheatmap", 
  "pROC", "RColorBrewer", "reshape2", "shiny", 
  "shinydashboard", "survival", "survminer", "tidyverse"
))
```

3. Launch the app:

```r
shiny::runApp("app.R")
```

---

## 🧪 Model Details

- **Algorithm:** LASSO-penalized Cox proportional hazards regression
- **Validation:** Cross-validation and external test sets (GSE28735, GSE62452)
- **Metric:** Concordance Index (C-index)

---

## 📊 Datasets Used

- [GSE28735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28735):  Contains paired pancreatic tumor and adjacent non-tumor tissue samples
- [GSE62452](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452):  Includes pancreatic ductal adenocarcinoma samples with survival information

The data is automatically downloaded from GEO during the first run of the application.
---

## 👥 Contributors

This application was developed by Biomed 8 team members:

- **Dylan**: Mathematics and Computer Science
- **Howard**: Data Science
- **Shuhan**: Data Science
- **Cherry**: Finance and Data Science
- **Ashley**: Data Science
- **Zhantao**: Data Science and Financial Mathematics

---
