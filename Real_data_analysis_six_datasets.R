
# ============================================================================
# REPRODUCIBILITY FIX (added 2026-08-08)
#
# real-data results to be reproducible and citable.
Sys.setenv(OMP_NUM_THREADS = "1")
options(rf.cores = 1, mc.cores = 1)
# ============================================================================

# ============================================================================
# Real_data_analysis_six_datasets.R
# Six real-life survival datasets:
# PBC, DBCD, CHOP, NKI70, VDX, TCGA Ovarian
# ============================================================================

rm(list = ls())

source("R/00_packages.R")
load_required_packages()
source("R/02_common_functions.R")
source("R/03_fit_agp_vi.R")
source("R/04b_fit_gl_blended_igh_vi.R")
source("R/05_fit_existing_methods.R")
source("R/06_evaluation.R")
source("R/08_plots.R")

dir.create("results/real_data", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/real_data", recursive = TRUE, showWarnings = FALSE)

old_results_files <- list.files("results/real_data", full.names = TRUE, recursive = TRUE)
if (length(old_results_files) > 0) {
  file.remove(old_results_files)
  cat("Cleared", length(old_results_files), "old file(s) from results/real_data/\n")
}
old_figure_files <- list.files("figures/real_data", full.names = TRUE, recursive = TRUE)
if (length(old_figure_files) > 0) {
  file.remove(old_figure_files)
  cat("Cleared", length(old_figure_files), "old file(s) from figures/real_data/\n")
}

extra_cran <- c(
  "survival", "randomForestSRC", "timeROC", "glmnet",
  "ggplot2", "dplyr", "tidyr", "penalized"
)

for (pkg in extra_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

extra_bioc <- c("Biobase", "breastCancerVDX", "curatedOvarianData")

for (pkg in extra_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  library(pkg, character.only = TRUE)
}

set.seed(2026)

train_fraction <- 0.70
rsf_quantile_cutoff <- 0.95
rsf_min_p <- 20
rsf_max_p <- 1000
rsf_ntree <- 500

absolute_threshold <- 0.15
z_threshold <- 2.33
selection_rule <- "or"
selection_method <- "auto"  # covariance methods: 95% CI; penalized Cox: native nonzero
credible_level <- 0.95

vi_max_iter <- 1000
gl_vi_max_iter <- 1000
gl_restart_on_failure <- TRUE
gl_stall_patience <- 30

gl_multistart_mode <- Sys.getenv("GL_MULTISTART_MODE", "adaptive")
if (!gl_multistart_mode %in% c("adaptive", "always")) {
  gl_multistart_mode <- "adaptive"
}
gl_temper_scales <- c(0, 0.25, 0.50, 0.75, 1.00)
vi_block_size <- 10

bayes_chains <- 2
bayes_iter <- 500
bayes_warmup <- 250

safe_dataset_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

make_numeric_matrix <- function(X) {
  X <- as.data.frame(X)
  X <- X[, vapply(X, is.numeric, logical(1)), drop = FALSE]
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  keep <- apply(X, 2, function(z) {
    zf <- z[is.finite(z)]
    length(zf) >= 2 && is.finite(stats::sd(zf)) && stats::sd(zf) > 0
  })
  X[, keep, drop = FALSE]
}

fit_train_preprocessor <- function(X_train) {
  X_train <- make_numeric_matrix(X_train)
  med <- apply(X_train, 2, function(z) stats::median(z[is.finite(z)], na.rm = TRUE))
  X_imp <- X_train
  for (j in seq_len(ncol(X_imp))) {
    bad <- !is.finite(X_imp[, j])
    if (any(bad)) X_imp[bad, j] <- med[j]
  }
  mu <- colMeans(X_imp)
  sig <- apply(X_imp, 2, stats::sd)
  keep <- is.finite(sig) & sig > 0 & is.finite(mu) & is.finite(med)
  list(columns = colnames(X_imp)[keep], median = med[keep], mean = mu[keep], sd = sig[keep])
}

apply_train_preprocessor <- function(X, prep) {
  X <- as.matrix(X[, prep$columns, drop = FALSE])
  storage.mode(X) <- "double"
  for (j in seq_len(ncol(X))) {
    bad <- !is.finite(X[, j])
    if (any(bad)) X[bad, j] <- prep$median[j]
  }
  X <- sweep(X, 2, prep$mean, "-")
  X <- sweep(X, 2, prep$sd, "/")
  X[!is.finite(X)] <- 0
  X
}

make_split <- function(X, time, status, train_fraction = 0.70, seed = 2026) {
  set.seed(seed)
  n <- nrow(X)
  train_id <- sample(seq_len(n), size = floor(train_fraction * n))

  test_id <- setdiff(seq_len(n), train_id)
  list(
    train_id = train_id,
    test_id = test_id,
    X_train = X[train_id, , drop = FALSE],
    X_test = X[test_id, , drop = FALSE],
    time_train = time[train_id],
    time_test = time[test_id],
    status_train = status[train_id],
    status_test = status[test_id]
  )
}

risk_score <- function(fit, X) {
  beta <- fit$beta
  if (is.null(beta)) return(rep(NA_real_, nrow(X)))
  beta[!is.finite(beta)] <- 0
  as.numeric(as.matrix(X) %*% beta)
}

screen_by_rsf_quantile <- function(X,
                                   time,
                                   status,
                                   quantile_cutoff = 0.95,
                                   min_p = 20,
                                   max_p = 1000,
                                   ntree = 500,
                                   seed = 2026) {
  set.seed(seed)

  X <- make_numeric_matrix(X)
  original_names <- colnames(X)
  if (is.null(original_names)) original_names <- paste0("V", seq_len(ncol(X)))
  safe_names <- paste0("X", seq_len(ncol(X)))
  colnames(X) <- safe_names

  name_map <- data.frame(
    Internal_Name = safe_names,
    Original_Variable = original_names,
    stringsAsFactors = FALSE
  )

  T <- as.numeric(time)
  event <- as.integer(status)
  train_data_rf <- data.frame(T = T, event = event, X, check.names = FALSE)

  # randomForestSRC parses survival formulas internally and expects the
  # response call to be exactly Surv(...), not survival::Surv(...).  The
  # namespace-qualified form can trigger parseFormula() errors on current
  # randomForestSRC versions.
  rf_model <- randomForestSRC::rfsrc(
    Surv(T, event) ~ .,
    data = train_data_rf,
    ntree = ntree,
    importance = TRUE,
    na.action = "na.impute",
    seed = -abs(as.integer(seed))
  )

  importance_scores <- randomForestSRC::vimp(rf_model)$importance
  importance_scores[!is.finite(importance_scores)] <- 0
  importance_scores <- importance_scores[names(importance_scores) %in% safe_names]

  cutoff_importance <- stats::quantile(
    importance_scores,
    quantile_cutoff,
    na.rm = TRUE,
    names = FALSE
  )

  selected_internal <- names(importance_scores[importance_scores > cutoff_importance])
  ranked_internal <- names(sort(importance_scores, decreasing = TRUE))

  if (length(selected_internal) < min_p) {
    selected_internal <- ranked_internal[seq_len(min(min_p, length(ranked_internal)))]
  }
  if (length(selected_internal) > max_p) {
    selected_internal <- ranked_internal[seq_len(min(max_p, length(ranked_internal)))]
  }
  selected_internal <- intersect(selected_internal, safe_names)
  if (length(selected_internal) == 0) {
    selected_internal <- ranked_internal[seq_len(min(min_p, length(ranked_internal)))]
  }

  selected_original <- name_map$Original_Variable[
    match(selected_internal, name_map$Internal_Name)
  ]

  cat("\nRSF screening summary\n")
  cat("Original p:", ncol(X), "\n")
  cat("Selected p:", length(selected_original), "\n")
  cat("Quantile cutoff:", quantile_cutoff, "\n")
  cat("Importance cutoff:", as.numeric(cutoff_importance), "\n")

  list(
    selected_variables = selected_original,
    selected_internal = selected_internal,
    importance = unname(importance_scores[selected_internal]),
    all_importance = importance_scores,
    name_map = name_map,
    threshold_rule = paste0("Training-only RSF importance > ", quantile_cutoff, " quantile"),
    threshold_value = as.numeric(cutoff_importance),
    original_p = ncol(X),
    selected_p = length(selected_original),
    rsf_fit = rf_model
  )
}

calc_time_auc <- function(method_name, fit, X_test, time_test, status_test, times) {
  lp <- risk_score(fit, X_test)

  out <- tryCatch(
    timeROC::timeROC(
      T = time_test,
      delta = status_test,
      marker = lp,
      cause = 1,
      weighting = "marginal",
      times = times,
      iid = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(data.frame(Method = method_name, Time = times, AUC = NA_real_))
  }

  data.frame(Method = method_name, Time = times, AUC = as.numeric(out$AUC))
}

evaluate_real_model <- function(fit,
                                X_train,
                                time_train,
                                status_train,
                                X_test,
                                time_test,
                                status_test,
                                eval_time = NULL,
                                threshold = 0.15,
                                z_threshold = 2.33,
                                selection_rule = "or",
                                selection_method = "auto",
                                credible_level = 0.95,
                                auc_times = NULL) {
  beta_hat <- fit$beta

  if (is.null(beta_hat) || length(beta_hat) != ncol(X_test) || all(is.na(beta_hat))) {
    return(
      data.frame(
        Method = fit$method,
        C_index = NA_real_,
        Brier_Score = NA_real_,
        Mean_Time_Dependent_AUC = NA_real_,
        Mean_Selected = NA_real_
      )
    )
  }

  beta_hat[!is.finite(beta_hat)] <- 0
  risk <- predict_risk(X = X_test, beta = beta_hat)

  if (is.null(eval_time) || !is.finite(eval_time)) {
    eval_time <- if (sum(status_train == 1) > 0) {
      as.numeric(stats::quantile(time_train[status_train == 1], 0.50, na.rm = TRUE))
    } else {
      as.numeric(stats::median(time_train, na.rm = TRUE))
    }
  }

  basehaz <- baseline_hazard_breslow(
    X = X_train,
    time = time_train,
    status = status_train,
    beta = beta_hat
  )

  surv <- predict_survival(
    X = X_test,
    beta = beta_hat,
    basehaz = basehaz,
    eval_time = eval_time
  )

  selection <- select_variables(
    fit = fit,
    absolute_threshold = threshold,
    z_threshold = z_threshold,
    rule = selection_rule,
    selection_method = selection_method,
    credible_level = credible_level
  )

  if (is.null(auc_times)) {
    auc_times <- as.numeric(
      stats::quantile(
        time_test[status_test == 1],
        probs = c(0.25, 0.50, 0.75),
        na.rm = TRUE
      )
    )
    auc_times <- unique(auc_times[is.finite(auc_times)])
  }

  if (length(auc_times) < 2) {
    auc_times <- as.numeric(
      stats::quantile(
        time_test,
        probs = c(0.25, 0.50, 0.75),
        na.rm = TRUE
      )
    )
    auc_times <- unique(auc_times[is.finite(auc_times)])
  }

  auc_df <- calc_time_auc(
    method_name = fit$method,
    fit = fit,
    X_test = X_test,
    time_test = time_test,
    status_test = status_test,
    times = auc_times
  )

  data.frame(
    Method = fit$method,
    C_index = harrell_cindex(time = time_test, status = status_test, risk = risk),
    Brier_Score = ipcw_brier(
      time = time_test,
      status = status_test,
      surv_prob = surv,
      eval_time = eval_time
    ),
    Mean_Time_Dependent_AUC = mean(auc_df$AUC, na.rm = TRUE),
    Mean_Selected = as.numeric(selection$n_selected)
  )
}

plot_auc_panel <- function(auc_results, datasets, save_path, title_text) {
  plot_data <- auc_results |>
    dplyr::filter(Dataset %in% datasets)

  plot_data$Dataset <- factor(plot_data$Dataset, levels = datasets)

  plot_data$Method <- factor(
    plot_data$Method,
    levels = c(
      "Cox-Lasso",
      "Cox-ElasticNet",
      "Bayesian-Lasso-Cox",
      "AGP-VI",
      "Global-Local Gaussian–IGH-VI"
    )
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Time,
      y = AUC,
      color = Method,
      linetype = Method,
      group = Method
    )
  ) +
    ggplot2::geom_line(linewidth = 1.15, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.4, na.rm = TRUE) +
    ggplot2::facet_wrap(~ Dataset, nrow = 1, scales = "free_x") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(
      title = title_text,
      x = "Evaluation time",
      y = "Time-dependent AUC",
      color = "Method",
      linetype = "Method"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
      strip.text = ggplot2::element_text(face = "bold", size = 13),
      axis.title = ggplot2::element_text(face = "bold", size = 13),
      axis.text = ggplot2::element_text(face = "bold", size = 11),
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(save_path, p, width = 14, height = 5.5, dpi = 600)

  ggplot2::ggsave(
    sub("\\.png$", ".pdf", save_path),
    p,
    width = 14,
    height = 5.5,
    device = cairo_pdf
  )

  p
}

load_pbc_data <- function() {
  data(pbc, package = "survival")
  # Preserve the same complete-case PBC cohort used in v5/v9 (n = 276) so
  # results remain directly comparable across package versions.  This filter
  # is performed before outcome modelling and no full-data standardization is
  # applied; scaling is still learned from the training partition only.
  pbc <- na.omit(pbc)
  time <- as.numeric(pbc$time)
  status <- as.integer(pbc$status == 2)
  X <- pbc[, -c(1, 2, 3), drop = FALSE]
  X <- make_numeric_matrix(X)

  list(name = "PBC", X = X, time = time, status = status)
}

load_dbcd_data <- function(file_path = "DBCD_Data.csv") {
  if (!file.exists(file_path)) {
    stop("Cannot find ", file_path, ". Put DBCD_Data.csv in the project working directory.")
  }

  DBCD_Data <- read.csv(file_path)
  DBCD_Data$status <- as.integer(DBCD_Data$status)
  DBCD_Data$time <- as.numeric(DBCD_Data$time)

  X <- DBCD_Data[, -c(1, 2), drop = FALSE]
  X <- make_numeric_matrix(X)

  list(name = "DBCD", X = X, time = DBCD_Data$time, status = DBCD_Data$status)
}

load_chop_data <- function(file_path = "chop.csv") {
  if (!file.exists(file_path)) {
    stop("Cannot find ", file_path, ". Put chop.csv in the project working directory.")
  }

  Datachop <- read.csv(file_path)
  Datachop <- Datachop[, -1, drop = FALSE]
  Datachop$status <- as.integer(Datachop$status)
  Datachop$survtime <- as.numeric(Datachop$survtime)

  keep <- is.finite(Datachop$survtime) &
    Datachop$survtime > 0 &
    is.finite(Datachop$status)

  Datachop <- Datachop[keep, , drop = FALSE]

  X <- Datachop[, -c(1, 2), drop = FALSE]
  X <- make_numeric_matrix(X)

  list(name = "CHOP", X = X, time = Datachop$survtime, status = Datachop$status)
}

load_nki70 <- function() {
  data("nki70", package = "penalized")
  time <- as.numeric(nki70$time)
  status <- as.integer(nki70$event)
  X <- nki70[, !(names(nki70) %in% c("time", "event")), drop = FALSE]
  X <- make_numeric_matrix(X)

  list(name = "NKI70", X = X, time = time, status = status)
}

load_vdx <- function() {
  data("vdx", package = "breastCancerVDX")
  X <- t(Biobase::exprs(vdx))
  ph <- Biobase::pData(vdx)

  time_var <- grep("t.dmfs|time|survival", names(ph), ignore.case = TRUE, value = TRUE)[1]
  event_var <- grep("e.dmfs|event|status", names(ph), ignore.case = TRUE, value = TRUE)[1]

  time <- as.numeric(ph[[time_var]])
  status <- as.integer(ph[[event_var]])

  keep <- is.finite(time) &
    is.finite(status) &
    time > 0

  X <- X[keep, , drop = FALSE]
  time <- time[keep]
  status <- status[keep]
  X <- make_numeric_matrix(X)

  list(name = "VDX", X = X, time = time, status = status)
}


load_ovarian <- function() {
  data("TCGA_eset", package = "curatedOvarianData")
  
  eset <- TCGA_eset
  X <- t(Biobase::exprs(eset))
  ph <- Biobase::pData(eset)
  
  # Use curated overall survival variables if available
  if ("days_to_death" %in% names(ph)) {
    time <- as.numeric(ph$days_to_death)
  } else if ("os" %in% names(ph)) {
    time <- as.numeric(ph$os)
  } else {
    time_var <- grep(
      "days_to_death|os|survival|time",
      names(ph),
      ignore.case = TRUE,
      value = TRUE
    )[1]
    time <- as.numeric(ph[[time_var]])
  }
  
  if ("vital_status" %in% names(ph)) {
    raw_status <- ph$vital_status
  } else if ("vital_status_at_death" %in% names(ph)) {
    raw_status <- ph$vital_status_at_death
  } else if ("deceased" %in% names(ph)) {
    raw_status <- ph$deceased
  } else {
    status_var <- grep(
      "vital_status|deceased|death|event|status",
      names(ph),
      ignore.case = TRUE,
      value = TRUE
    )[1]
    raw_status <- ph[[status_var]]
  }
  
  raw_status_chr <- tolower(trimws(as.character(raw_status)))
  
  status <- ifelse(
    raw_status_chr %in% c(
      "dead",
      "deceased",
      "death",
      "1",
      "yes",
      "true"
    ),
    1,
    0
  )
  
  # If still no deaths, try common curatedOvarianData field names
  if (sum(status == 1, na.rm = TRUE) == 0) {
    possible_status <- grep(
      "os.*event|event.*os|death|deceased|vital",
      names(ph),
      ignore.case = TRUE,
      value = TRUE
    )
    
    for (sv in possible_status) {
      tmp <- ph[[sv]]
      
      if (is.numeric(tmp) || is.integer(tmp)) {
        tmp_status <- as.integer(tmp)
      } else {
        tmp_chr <- tolower(trimws(as.character(tmp)))
        tmp_status <- ifelse(
          tmp_chr %in% c(
            "dead",
            "deceased",
            "death",
            "1",
            "yes",
            "true"
          ),
          1,
          0
        )
      }
      
      if (sum(tmp_status == 1, na.rm = TRUE) > 0) {
        status <- tmp_status
        cat("Using TCGA ovarian status variable:", sv, "\n")
        break
      }
    }
  }
  
  keep <- is.finite(time) &
    is.finite(status) &
    time > 0
  
  X <- X[keep, , drop = FALSE]
  time <- time[keep]
  status <- status[keep]
  
  if (sum(status == 1, na.rm = TRUE) == 0) {
    stop(
      "TCGA ovarian still has no deaths after status recoding. ",
      "Run names(Biobase::pData(TCGA_eset)) to inspect survival variables."
    )
  }
  
  X <- make_numeric_matrix(X)
  
  cat("\n====================================\n")
  cat("Dataset: TCGA Ovarian\n")
  cat("Sample size (n):", nrow(X), "\n")
  cat("Number of predictors (p):", ncol(X), "\n")
  cat("Deaths/events:", sum(status == 1), "\n")
  cat("Censored:", sum(status == 0), "\n")
  cat("====================================\n")
  
  list(
    name = "TCGA Ovarian",
    X = X,
    time = time,
    status = status
  )
}

real_datasets <- list(
  load_pbc_data()
  load_dbcd_data("DBCD_Data.csv"),
  load_chop_data("chop.csv"),
  load_nki70(),
  load_vdx(),
  load_ovarian()
)

real_results <- list()
auc_results <- list()
counter <- 1


# ============================================================================
# Persist all fitted VI models from all six datasets
# ============================================================================

all_agp_fits <- list()
all_gl_fits <- list()

fit_dir <- "results/real_data/fits"
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)


for (dat in real_datasets ) {
  cat("\n========================================\n")
  cat("Real dataset:", dat$name, "\n")
  cat("Original n =", nrow(dat$X), " original p =", ncol(dat$X), "\n")
  cat("========================================\n")

  original_n <- nrow(dat$X)
  original_p <- ncol(dat$X)

  dat_keep <- is.finite(dat$time) &
    dat$time > 0 &
    is.finite(dat$status)

  dat$X <- dat$X[dat_keep, , drop = FALSE]
  dat$time <- dat$time[dat_keep]
  dat$status <- dat$status[dat_keep]

  cat("After removing invalid survival times:\n")
  cat("Remaining n =", nrow(dat$X), "\n")

  # Split FIRST. All supervised screening, imputation and scaling below use
  # the training set only so the test outcomes remain untouched.
  split_raw <- make_split(
    X = dat$X,
    time = dat$time,
    status = dat$status,
    train_fraction = train_fraction,
    seed = 2026
  )

  if (ncol(split_raw$X_train) <= 100) {
    selected_variables <- colnames(split_raw$X_train)
    rsf_screen <- list(
      selected_variables = selected_variables,
      importance = rep(NA_real_, length(selected_variables)),
      threshold_rule = "No screening: p <= 100",
      threshold_value = NA_real_,
      original_p = ncol(split_raw$X_train),
      selected_p = length(selected_variables)
    )
  } else {
    rsf_screen <- screen_by_rsf_quantile(
      X = split_raw$X_train,
      time = split_raw$time_train,
      status = split_raw$status_train,
      quantile_cutoff = rsf_quantile_cutoff,
      min_p = rsf_min_p,
      max_p = rsf_max_p,
      ntree = rsf_ntree,
      seed = 2026
    )
    selected_variables <- rsf_screen$selected_variables
  }

  # Apply the training-selected variables to both partitions, then fit a
  # median-imputation + standardization recipe on training predictors only.
  X_train_selected <- split_raw$X_train[, selected_variables, drop = FALSE]
  X_test_selected <- split_raw$X_test[, selected_variables, drop = FALSE]
  prep <- fit_train_preprocessor(X_train_selected)

  split <- split_raw
  split$X_train <- apply_train_preprocessor(X_train_selected, prep)
  split$X_test <- apply_train_preprocessor(X_test_selected, prep)

  # If preprocessing removed a zero-variance training variable, keep the RSF
  # metadata aligned with the columns actually fitted by the models.
  selected_variables_final <- colnames(split$X_train)
  keep_rsf <- match(selected_variables_final, rsf_screen$selected_variables)
  selected_importance_final <- rsf_screen$importance[keep_rsf]

  cat("\nAfter training-only RSF/preprocessing:\n")
  cat("Selected p =", ncol(split$X_train), "\n")
  cat("Threshold rule =", rsf_screen$threshold_rule, "\n")
  cat("Threshold value =", rsf_screen$threshold_value, "\n")

  write.csv(
    data.frame(
      Variable = selected_variables_final,
      Importance = as.numeric(selected_importance_final)
    ),
    paste0("results/real_data/RSF_selected_variables_", safe_dataset_name(dat$name), ".csv"),
    row.names = FALSE
  )

  # Evaluation horizons are determined from TRAINING survival times, never
  # from test outcomes.
  times_auc <- as.numeric(
    quantile(
      split$time_train[split$status_train == 1],
      probs = c(0.25, 0.50, 0.75),
      na.rm = TRUE
    )
  )
  times_auc <- unique(times_auc[is.finite(times_auc)])

  if (length(times_auc) < 2) {
    times_auc <- as.numeric(
      quantile(
        split$time_train,
        probs = c(0.25, 0.50, 0.75),
        na.rm = TRUE
      )
    )
  }

  t0 <- proc.time()[3]
  fit_lasso <- fit_cox_lasso(split$X_train, split$time_train, split$status_train)
  time_lasso <- proc.time()[3] - t0

  t0 <- proc.time()[3]
  fit_enet <- fit_cox_elasticnet(split$X_train, split$time_train, split$status_train, alpha = 0.5)
  time_enet <- proc.time()[3] - t0

  t0 <- proc.time()[3]
  fit_bayes <- fit_bayesian_lasso_cox(
    split$X_train,
    split$time_train,
    split$status_train,
    chains = bayes_chains,
    iter = bayes_iter,
    warmup = bayes_warmup,
    seed = 2026
  )
  time_bayes <- proc.time()[3] - t0

  t0 <- proc.time()[3]
  fit_agp <- fit_agp_vi(
    split$X_train,
    split$time_train,
    split$status_train,
    max_iter = vi_max_iter,
    eta = 0.001,
    block_size = vi_block_size
  )
  time_agp <- proc.time()[3] - t0

  t0 <- proc.time()[3]

  # general adaptive  GL-IGH-VI.
  fit_gl <- fit_gl_blended_igh_vi_adaptive(
    X = split$X_train,
    time = split$time_train,
    status = split$status_train,
    agp_fit = fit_agp,
    temper_scales = gl_temper_scales,
    multistart_mode = gl_multistart_mode,
    fallback_fraction_trigger = 0.10,
    verbose_multistart = TRUE,
    max_iter = gl_vi_max_iter,
    eta = 1.0,
    tol = 1e-4,
    patience = 5,
    stall_patience = gl_stall_patience,
    block_size = vi_block_size,
    rescue_mean_only = TRUE,
    rescue_raw_gradient = TRUE
  )

  time_gl <- proc.time()[3] - t0

  # Save initialization diagnostics for every dataset.
  if (!is.null(fit_gl$multistart_diagnostics)) {
    write.csv(
      fit_gl$multistart_diagnostics,
      paste0(
        "results/real_data/GL_IGH_multistart_",
        safe_dataset_name(dat$name),
        ".csv"
      ),
      row.names = FALSE
    )
  }

  # Posterior coefficient diagnostics: useful for detecting artificial
  # plateaus and for auditing credible-interval selection.
  gl_se <- sqrt(pmax(diag(fit_gl$covariance), 1e-12))
  gl_z <- abs(fit_gl$beta) / gl_se
  gl_zcrit <- stats::qnorm(1 - (1 - credible_level) / 2)
  gl_diag <- data.frame(
    Variable = colnames(split$X_train),
    Beta = as.numeric(fit_gl$beta),
    Posterior_SE = gl_se,
    CI_Lower = fit_gl$beta - gl_zcrit * gl_se,
    CI_Upper = fit_gl$beta + gl_zcrit * gl_se,
    Z_Abs = gl_z,
    Selected_95CI = (fit_gl$beta - gl_zcrit * gl_se > 0) | (fit_gl$beta + gl_zcrit * gl_se < 0),
    stringsAsFactors = FALSE
  )
  write.csv(gl_diag,
            paste0("results/real_data/GL_IGH_posterior_diagnostics_", safe_dataset_name(dat$name), ".csv"),
            row.names = FALSE)

  # DIAGNOSTIC (added to investigate the Mean_Selected=0 anomaly seen on
  # CHOP/DBCD/NKI70
  
  cat(sprintf(
    paste0(
      "  [diagnostic] %s GL-IGH-VI: status=%s, p=%d, iterations=%d, ",
      "n_fallback_used=%d, max_failed_streak=%d, multistart=%s, selected_start=%s
"
    ),
    dat$name,
    fit_gl$convergence_status,
    ncol(split$X_train),
    fit_gl$iterations,
    fit_gl$n_fallback_used,
    fit_gl$max_failed_step_streak,
    ifelse(isTRUE(fit_gl$multistart_used), "TRUE", "FALSE"),
    ifelse(
      is.null(fit_gl$selected_start_label),
      "cold",
      fit_gl$selected_start_label
    )
  ))

  print_elbo_diagnostics(fit_gl, "GL-IGH-VI", dat$name)


  # single-panel plot.
  plot_vi_elbo_two_panel(
    fit_gl,
    paste(dat$name, "GL-IGH-VI"),
    paste0("figures/real_data/", safe_dataset_name(dat$name), "_GL_IGH_ELBO_two_panel.png")
  )

  fit_lasso$method <- "Cox-Lasso"
  fit_enet$method <- "Cox-ElasticNet"
  fit_bayes$method <- "Bayesian-Lasso-Cox"
  fit_agp$method <- "AGP-VI"
  fit_gl$method <- "Global-Local Gaussian–IGH-VI"


  # ------------------------------------------------------------------------
  # Preserve current dataset fits before the loop advances
  # ------------------------------------------------------------------------

  dataset_key <- safe_dataset_name(dat$name)

  all_agp_fits[[dat$name]] <- fit_agp
  all_gl_fits[[dat$name]] <- fit_gl

  saveRDS(
    fit_agp,
    file.path(fit_dir, paste0(dataset_key, "_AGP_VI_fit.rds"))
  )

  saveRDS(
    fit_gl,
    file.path(fit_dir, paste0(dataset_key, "_GL_IGH_VI_fit.rds"))
  )


  fits <- list(
    "Cox-Lasso" = fit_lasso,
    "Cox-ElasticNet" = fit_enet,
    "Bayesian-Lasso-Cox" = fit_bayes,
    "AGP-VI" = fit_agp,
    "Global-Local Gaussian–IGH-VI" = fit_gl
  )

  runtimes <- c(
    "Cox-Lasso" = time_lasso,
    "Cox-ElasticNet" = time_enet,
    "Bayesian-Lasso-Cox" = time_bayes,
    "AGP-VI" = time_agp,
    "Global-Local Gaussian–IGH-VI" = time_gl
  )

  for (method_name in names(fits)) {
    result <- evaluate_real_model(
      fit = fits[[method_name]],
      X_train = split$X_train,
      time_train = split$time_train,
      status_train = split$status_train,
      X_test = split$X_test,
      time_test = split$time_test,
      status_test = split$status_test,
      threshold = absolute_threshold,
      z_threshold = z_threshold,
      selection_rule = selection_rule,
      selection_method = selection_method,
      credible_level = credible_level,
      auc_times = times_auc
    )

    auc_df <- calc_time_auc(
      method_name = method_name,
      fit = fits[[method_name]],
      X_test = split$X_test,
      time_test = split$time_test,
      status_test = split$status_test,
      times = times_auc
    )

    auc_df$Dataset <- dat$name
    auc_results[[counter]] <- auc_df

    result$Dataset <- dat$name
    result$Original_n <- original_n
    result$Original_p <- original_p
    result$RSF_Selected_p <- ncol(split$X_train)
    result$RSF_Threshold_Rule <- rsf_screen$threshold_rule
    result$RSF_Threshold_Value <- as.numeric(rsf_screen$threshold_value)
    result$Time_s <- as.numeric(unname(runtimes[method_name]))
    result$Convergence_Status <- if (!is.null(fits[[method_name]]$convergence_status)) fits[[method_name]]$convergence_status else NA_character_
    result$Fallbacks <- if (!is.null(fits[[method_name]]$n_fallback_used)) fits[[method_name]]$n_fallback_used else NA_integer_
    result$Max_Failed_Step_Streak <- if (!is.null(fits[[method_name]]$max_failed_step_streak)) fits[[method_name]]$max_failed_step_streak else NA_integer_
    result$Adaptive_Restart <- if (identical(method_name, "Global-Local Gaussian–IGH-VI")) {
      if (!is.null(fits[[method_name]]$adaptive_restart_used)) {
        isTRUE(fits[[method_name]]$adaptive_restart_used)
      } else {
        FALSE
      }
    } else {
      NA
    }

    result$Multistart_Used <- if (identical(method_name, "Global-Local Gaussian–IGH-VI")) {
      if (!is.null(fits[[method_name]]$multistart_used)) {
        isTRUE(fits[[method_name]]$multistart_used)
      } else {
        FALSE
      }
    } else {
      NA
    }

    result$Selected_Start <- if (identical(method_name, "Global-Local Gaussian–IGH-VI")) {
      if (!is.null(fits[[method_name]]$selected_start_label)) {
        as.character(fits[[method_name]]$selected_start_label)
      } else {
        "cold"
      }
    } else {
      NA_character_
    }
    sel_tmp <- select_variables(fits[[method_name]], selection_method = selection_method, credible_level = credible_level)
    result$Selection_Method <- sel_tmp$selection_method
    result$Credible_Level <- sel_tmp$credible_level

    result <- result[
      ,
      c(
        "Dataset",
        "Original_n",
        "Original_p",
        "RSF_Selected_p",
        "RSF_Threshold_Rule",
        "RSF_Threshold_Value",
        "Method",
        "C_index",
        "Brier_Score",
        "Mean_Time_Dependent_AUC",
        "Mean_Selected",
        "Time_s",
        "Convergence_Status",
        "Fallbacks",
        "Max_Failed_Step_Streak",
        "Adaptive_Restart",
        "Multistart_Used",
        "Selected_Start",
        "Selection_Method",
        "Credible_Level"
      )
    ]

    real_results[[counter]] <- result
    counter <- counter + 1
  }

  cat(
    "\n[PLOT CHECK]", dat$name,
    "AGP-VI:",
    "iterations =", length(fit_agp$elbo),
    "final ELBO =", tail(fit_agp$elbo, 1),
    "status =",
    ifelse(is.null(fit_agp$convergence_status), "unknown", fit_agp$convergence_status),
    "\n"
  )

  cat(
    "[PLOT CHECK]", dat$name,
    "GL-IGH-VI:",
    "iterations =", length(fit_gl$elbo),
    "final ELBO =", tail(fit_gl$elbo, 1),
    "status =", fit_gl$convergence_status,
    "start =",
    ifelse(is.null(fit_gl$selected_start_label), "cold", fit_gl$selected_start_label),
    "\n"
  )

  plot_vi_elbo(
    fit_agp,
    paste(dat$name, "AGP-VI"),
    paste0("figures/real_data/", safe_dataset_name(dat$name), "_AGP_ELBO.png")
  )

  plot_vi_elbo(
    fit_gl,
    paste(dat$name, "GL-IGH-VI"),
    paste0("figures/real_data/", safe_dataset_name(dat$name), "_GL_IGH_ELBO.png")
  )

  gc()
}


# Save collections of all six VI fits for later inspection
saveRDS(
  all_agp_fits,
  file.path(fit_dir, "ALL_AGP_VI_fits.rds")
)

saveRDS(
  all_gl_fits,
  file.path(fit_dir, "ALL_GL_IGH_VI_fits.rds")
)


real_results <- dplyr::bind_rows(real_results)
auc_results <- do.call(rbind, auc_results)

real_summary_results <- real_results |>
  dplyr::group_by(Dataset, Method) |>
  dplyr::summarise(
    C_index = mean(C_index, na.rm = TRUE),
    Brier_Score = mean(Brier_Score, na.rm = TRUE),
    Mean_Time_Dependent_AUC = mean(Mean_Time_Dependent_AUC, na.rm = TRUE),
    Mean_Selected = mean(Mean_Selected, na.rm = TRUE),
    Time_s = mean(Time_s, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(real_results, "results/real_data/real_data_raw_results.csv", row.names = FALSE)
write.csv(real_summary_results, "results/real_data/real_data_summary_results.csv", row.names = FALSE)
write.csv(auc_results, "results/real_data/time_dependent_auc_results.csv", row.names = FALSE)

latex_real_table <- real_summary_results |>
  dplyr::mutate(
    C_index_num = round(C_index, 4),
    Brier_Score_num = round(Brier_Score, 4),
    AUC_num = round(Mean_Time_Dependent_AUC, 4),
    Mean_Selected_num = round(Mean_Selected, 0)
  ) |>
  dplyr::group_by(Dataset) |>
  dplyr::mutate(
    Best_C_index = is.finite(C_index) & abs(C_index - max(C_index, na.rm = TRUE)) < 1e-12,
    Best_Brier = is.finite(Brier_Score) & abs(Brier_Score - min(Brier_Score, na.rm = TRUE)) < 1e-12,
    Best_AUC = is.finite(Mean_Time_Dependent_AUC) & abs(Mean_Time_Dependent_AUC - max(Mean_Time_Dependent_AUC, na.rm = TRUE)) < 1e-12
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    C_index = ifelse(Best_C_index, paste0("\\textbf{", sprintf("%.4f", C_index_num), "}"), sprintf("%.4f", C_index_num)),
    Brier_Score = ifelse(Best_Brier, paste0("\\textbf{", sprintf("%.4f", Brier_Score_num), "}"), sprintf("%.4f", Brier_Score_num)),
    Mean_Time_Dependent_AUC = ifelse(Best_AUC, paste0("\\textbf{", sprintf("%.4f", AUC_num), "}"), sprintf("%.4f", AUC_num)),
    Mean_Selected = as.character(Mean_Selected_num)
  )

latex_lines <- c(
  "\\renewcommand{\\arraystretch}{1.4}",
  "\\begin{table}[ht]",
  "\\centering",
  "\\small",
  "\\caption{Real-life survival analysis results across six datasets.}",
  "\\resizebox{\\textwidth}{!}{",
  "\\begin{tabular}{|l|l|c|c|c|c|}",
  "\\hline",
  "Dataset & Method & C-index & Brier Score & Mean Time-Dependent AUC & No. Selected \\\\",
  "\\hline"
)

for (i in seq_len(nrow(latex_real_table))) {
  row_i <- latex_real_table[i, ]
  latex_lines <- c(
    latex_lines,
    paste(
      row_i$Dataset,
      "&",
      row_i$Method,
      "&",
      row_i$C_index,
      "&",
      row_i$Brier_Score,
      "&",
      row_i$Mean_Time_Dependent_AUC,
      "&",
      row_i$Mean_Selected,
      "\\\\"
    )
  )

  if (i == nrow(latex_real_table) ||
      latex_real_table$Dataset[i] != latex_real_table$Dataset[i + 1]) {
    latex_lines <- c(latex_lines, "\\hline")
  }
}

latex_lines <- c(
  latex_lines,
  "\\end{tabular}",
  "}",
  "\\label{tab:real_data_results}",
  "\\end{table}"
)

writeLines(latex_lines, con = "results/real_data/real_data_summary_table.tex")

plot_auc_panel(
  auc_results = auc_results,
  datasets = c("PBC", "DBCD", "CHOP"),
  save_path = "figures/real_data/time_dependent_auc_panel_1_PBC_DBCD_CHOP.png",
  title_text = "Time-dependent AUC: PBC, DBCD, and CHOP"
)

plot_auc_panel(
  auc_results = auc_results,
  datasets = c("NKI70", "VDX", "TCGA Ovarian"),
  save_path = "figures/real_data/time_dependent_auc_panel_2_NKI70_VDX_TCGA.png",
  title_text = "Time-dependent AUC: NKI70, VDX, and TCGA Ovarian"
)

capture.output(sessionInfo(), file = "results/real_data/session_info.txt")

cat("\nReal-data analysis complete.\n")
cat("\nSummary CSV saved to: results/real_data/real_data_summary_results.csv\n")
cat("LaTeX table saved to: results/real_data/real_data_summary_table.tex\n")
cat("Figure 1 saved to: figures/real_data/time_dependent_auc_panel_1_PBC_DBCD_CHOP.png\n")
cat("Figure 2 saved to: figures/real_data/time_dependent_auc_panel_2_NKI70_VDX_TCGA.png\n\n")

print(real_summary_results)
