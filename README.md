#  A Hierarchical Bayesian Cox Model with Variational Inference for High-Dimensional Survival Analysis

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
├──run_simulation_four_methods
├──results/
├── figures/
├── 10_real_data_analysis_six_datasets.R
├── chop.csv
├── results/
├── figures/

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

# Quick Run

## 1. Clone the Repository

```bash
git clone https://github.com/Dydx1989/AGP-VI-and-Gaussian-IGH-VI.git
cd AGP-VI-and-Gaussian-IGH-VI
```

## 2. Install Required Packages

```r
source("R/00_packages.R")
load_required_packages()
```

## 3.Quick run

```r
source("R/00_packages.R")
load_required_packages()

source("R/01_simulate_data.R")
source("R/02_common_functions.R")
source("R/03_fit_agp_vi.R")
source("R/04b_fit_gl_blended_igh_vi.R")
source("R/05_fit_existing_methods.R")
source("R/06_evaluation.R")

set.seed(1)

dat <- simulate_cox_data(
  n = 100,
  p = 200,
  n_signal = 10,
  rho = 0.30,
  beta_value = 0.7,
  censoring_rate = 0.50,
  seed = 1
)
split <- train_test_split(
  dat$X,
  dat$time,
  dat$status,
  train_fraction = 0.70,
  seed = 1
)

eval_time <- as.numeric(
  quantile(
    split$time_train[split$status_train == 1],
    0.50
  )
)

fit_cox_lasso <- fit_cox_lasso(
  split$X_train,
  split$time_train,
  split$status_train
)

fit_elasticnet <- fit_cox_elasticnet(
  split$X_train,
  split$time_train,
  split$status_train,
  alpha = 0.5
)

fit_bayes_lasso <- fit_bayesian_lasso_cox(
  split$X_train,
  split$time_train,
  split$status_train,
  chains = 2,
  iter = 500,
  warmup = 250,
  seed = 1
)

fit_agp <- fit_agp_vi(
  split$X_train,
  split$time_train,
  split$status_train,
  max_iter = 50,
  block_size = 10
)

fit_gl_igh <- fit_gl_blended_igh_vi(
  split$X_train,
  split$time_train,
  split$status_train,
  max_iter = 50,
  block_size = 10
)

res_cox_lasso <- evaluate_model(
  fit_cox_lasso,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_elasticnet <- evaluate_model(
  fit_elasticnet,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_bayes_lasso <- evaluate_model(
  fit_bayes_lasso,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_agp <- evaluate_model(
  fit_agp,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_gl_igh <- evaluate_model(
  fit_gl_igh,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_cox_lasso$Method <- "Cox-Lasso"
res_elasticnet$Method <- "Cox-ElasticNet"
res_bayes_lasso$Method <- "Bayesian-Lasso-Cox"
res_agp$Method <- "AGP-VI"
res_gl_igh$Method <- "Global-Local Gaussian–IGH-VI"

quick_results <- rbind(
  res_cox_lasso,
  res_elasticnet,
  res_bayes_lasso,
  res_agp,
  res_gl_igh
)

print(quick_results)

```





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
