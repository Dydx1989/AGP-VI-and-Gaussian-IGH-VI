#  A Robust Hierarchical Bayesian Cox Model with Variational Inference for High-Dimensional Survival Analysis

Robust hierarchical Bayesian variational inference methods for high-dimensional Cox survival models, including AGP-VI and Global-Local Gaussian–IGH-VI with simulation studies, real genomic survival datasets, RSF feature screening, ELBO diagnostics, and time-dependent AUC analysis.

## Overview

This repository contains the complete reproducibility materials for the paper:

> **A Robust Hierarchical Bayesian Cox Model with Variational Inference for High-Dimensional Survival Analysis**

The project develops scalable hierarchical Bayesian variational inference methods for high-dimensional Cox proportional hazards models under sparse survival settings.

The repository includes:

* Full simulation studies
* Real-life high-dimensional survival data analyses
* Variational inference implementations
* Benchmark survival methods
* Time-dependent AUC and Integrated Brier Score evaluation
* ELBO convergence diagnostics
* Publication-quality figures and LaTeX tables

---

# Implemented Methods

The following methods are implemented and compared throughout the paper:

1. **Cox-Lasso**
2. **Cox-ElasticNet**
3. **Bayesian-Lasso-Cox**
4. **AGP-VI** (Adaptive Gaussian Prior Variational Inference)
5. **Global-Local Gaussian–IGH-VI**

---

# Repository Structure

```text
.
├── R/
│   ├── 00_packages.R
│   ├── 01_simulate_data.R
│   ├── 02_common_functions.R
│   ├── 03_fit_agp_vi.R
│   ├── 04_fit_igh_vi.R
│   ├── 04b_fit_gl_blended_igh_vi.R
│   ├── 05_fit_existing_methods.R
│   ├── 06_evaluation.R
│   ├── 07_summary_latex.R
│   └── 08_plots.R
│
├── results/
├── figures/
├── DBCD_Data.csv
├── chop.csv
├── 10_real_data_analysis_six_datasets.R
├── README.md
└── README_six_real_data_analysis.md
```

---

# Methodological Framework

The proposed framework introduces robust hierarchical Bayesian Cox proportional hazards models with scalable variational inference for high-dimensional survival analysis.

The methodology is designed to address:

* High-dimensional covariate spaces
* Sparse signal recovery
* Correlated predictors
* Heavy censoring
* Computational scalability

The proposed variational inference procedures provide efficient posterior approximation while maintaining strong predictive performance and variable selection accuracy.

---

# Simulation Study

The simulation study evaluates the methods under:

* High-dimensional sparse settings
* Correlated covariates
* Different predictor dimensions
* Moderate and heavy censoring
* Variable selection performance
* Prediction accuracy

## Simulation Settings

Typical configurations include:

```r
n = 100
p = c(50, 100, 200, 300, 400, 500)
n_signal = 10
rho = 0.15
beta_value = 1.20
censoring_rate = 0.50
```

## Performance Metrics

The following metrics are reported:

* Harrell's C-index
* Brier Score
* Time-dependent AUC
* True Positive Rate (TPR)
* False Positive Rate (FPR)
* Variable selection accuracy
* Computational runtime

---

# Real-Life Survival Datasets

The repository includes analyses on six real-life survival datasets.

| Dataset      | Subjects (n) | Predictors (p) | Description                                        |
| ------------ | -----------: | -------------: | -------------------------------------------------- |
| PBC          |          276 |             16 | Primary Biliary Cirrhosis survival dataset         |
| DBCD         |          295 |          4,919 | Dutch Breast Cancer gene expression dataset        |
| CHOP         |          181 |          3,833 | Diffuse Large B-cell Lymphoma dataset              |
| NKI70        |          144 |             71 | Netherlands Cancer Institute breast cancer dataset |
| VDX          |          344 |         22,283 | High-dimensional breast cancer microarray dataset  |
| TCGA Ovarian |          557 |         13,104 | TCGA ovarian cancer genomic survival dataset       |

---

# RSF Feature Screening

For ultra-high-dimensional datasets, the repository performs feature screening using:

## Random Survival Forests (RSF)

The screening procedure:

1. Fits RSF on the complete covariate space
2. Computes variable importance scores
3. Selects variables exceeding a quantile importance threshold
4. Passes selected covariates into the survival models

Typical configuration:

```r
quantile_cutoff = 0.95
min_p = 20
max_p = 1000
ntree = 500
```

---

# Running the Code

## 1. Install Required CRAN Packages

```r
install.packages(c(
  "survival",
  "glmnet",
  "MASS",
  "ggplot2",
  "dplyr",
  "tidyr",
  "gridExtra",
  "randomForestSRC",
  "timeROC",
  "penalized"
))
```

## 2. Install Required Bioconductor Packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "Biobase",
  "breastCancerVDX",
  "curatedOvarianData"
))
```

---

# Running Simulation Experiments

```r
source("run_simulation_four_methods.R")
```

Outputs include:

* Simulation summary tables
* LaTeX tables
* ELBO convergence plots
* Time-dependent AUC figures
* Integrated Brier Score plots

---

# Running Real-Life Data Analyses

```r
source("10_real_data_analysis_six_datasets.R")
```

Outputs are automatically saved in:

```text
results/real_data/
figures/real_data/
```

---

# Output Files

The repository automatically generates:

## Tables

* Simulation summary tables
* Real-data performance tables
* Publication-ready LaTeX tables

## Figures

* ELBO convergence diagnostics
* Coefficient recovery plots
* Time-dependent AUC panels
* Integrated Brier Score plots

## CSV Results

* Raw performance metrics
* Time-dependent AUC values
* RSF-selected variables

---

# Reproducibility

All experiments are fully reproducible.

Set the random seed before execution:

```r
set.seed(2026)
```

All scripts used to generate the manuscript results are included in this repository.

---

# Citation

If you use this repository, please cite:

```text
Author(s). A Robust Hierarchical Bayesian Cox Model with
Variational Inference for High-Dimensional Survival Analysis.
```

---

# References

## Survival Package

Therneau, T. M. (2024).
*A Package for Survival Analysis in R.*

## NKI70 Breast Cancer Dataset

van't Veer, L. J., et al. (2002).
Gene expression profiling predicts clinical outcome of breast cancer.
Nature.

## TCGA Ovarian Dataset

The Cancer Genome Atlas Research Network (2011).
Integrated genomic analyses of ovarian carcinoma.
Nature.

## DLBCL / CHOP Dataset

Shipp, M. A., et al. (2002).
Diffuse large B-cell lymphoma outcome prediction by gene-expression profiling and supervised machine learning.
Nature Medicine.

---

# Contact

For questions regarding the code or reproducibility materials, please contact the repository author.
