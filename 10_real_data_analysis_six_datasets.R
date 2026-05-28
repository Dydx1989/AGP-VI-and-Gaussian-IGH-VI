
# ============================================================================
# 10_real_data_analysis_six_datasets.R
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

vi_max_iter <- 150
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
    stats::sd(z, na.rm = TRUE) > 0 && all(is.finite(z))
  })

  X[, keep, drop = FALSE]
}

standardize_x <- function(X) {
  X <- make_numeric_matrix(X)
  X <- scale(X)
  X[!is.finite(X)] <- 0
  X
}

make_split <- function(X, time, status, train_fraction = 0.70, seed = 2026) {
  set.seed(seed)
  n <- nrow(X)
  train_id <- sample(seq_len(n), size = floor(train_fraction * n))

  list(
    X_train = X[train_id, , drop = FALSE],
    X_test = X[-train_id, , drop = FALSE],
    time_train = time[train_id],
    time_test = time[-train_id],
    status_train = status[train_id],
    status_test = status[-train_id]
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
  
  X <- scale(data.matrix(X))
  X[!is.finite(X)] <- 0
  
  # Important: force safe predictor names
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  
  T <- as.numeric(time)
  event <- as.integer(status)
  
  train_data_rf <- data.frame(
    T = T,
    event = event,
    X,
    check.names = FALSE
  )
  
  rf_model <- rfsrc(
    Surv(T, event) ~ .,
    data = train_data_rf,
    ntree = ntree,
    importance = TRUE,
    na.action = "na.impute"
  )
  
  importance_scores <- randomForestSRC::vimp(rf_model)$importance
  importance_scores[!is.finite(importance_scores)] <- 0
  
  # Keep only variables that still match X column names
  importance_scores <- importance_scores[
    names(importance_scores) %in% colnames(X)
  ]
  
  cutoff_importance <- stats::quantile(
    importance_scores,
    quantile_cutoff,
    na.rm = TRUE
  )
  
  selected_features <- names(
    importance_scores[importance_scores > cutoff_importance]
  )
  
  ranked_features <- names(
    sort(
      importance_scores,
      decreasing = TRUE
    )
  )
  
  if (length(selected_features) < min_p) {
    selected_features <- ranked_features[
      seq_len(
        min(
          min_p,
          length(ranked_features)
        )
      )
    ]
  }
  
  if (length(selected_features) > max_p) {
    selected_features <- ranked_features[
      seq_len(max_p)
    ]
  }
  
  # Final protection against subscript errors
  selected_features <- intersect(
    selected_features,
    colnames(X)
  )
  
  if (length(selected_features) == 0) {
    selected_features <- ranked_features[
      seq_len(
        min(
          min_p,
          length(ranked_features)
        )
      )
    ]
  }
  
  X_filtered <- X[
    ,
    selected_features,
    drop = FALSE
  ]
  
  cat("\nRSF screening summary\n")
  cat("Original p:", ncol(X), "\n")
  cat("Selected p:", ncol(X_filtered), "\n")
  cat("Quantile cutoff:", quantile_cutoff, "\n")
  cat("Importance cutoff:", as.numeric(cutoff_importance), "\n")
  
  list(
    X_selected = X_filtered,
    selected_variables = selected_features,
    importance = importance_scores[selected_features],
    all_importance = importance_scores,
    threshold_rule = paste0(
      "RSF importance > ",
      quantile_cutoff,
      " quantile"
    ),
    threshold_value = as.numeric(cutoff_importance),
    original_p = ncol(X),
    selected_p = ncol(X_filtered),
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
    rule = selection_rule
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
  pbc <- na.omit(pbc)
  time <- as.numeric(pbc$time)
  status <- as.integer(pbc$status == 2)
  X <- pbc[, -c(1, 2, 3), drop = FALSE]
  X <- standardize_x(X)

  list(name = "PBC", X = X, time = time, status = status)
}

load_dbcd_data <- function(file_path = "DBCD_Data.csv") {
  if (!file.exists(file_path)) {
    stop("Cannot find ", file_path, ". Put DBCD_Data.csv in the project working directory.")
  }

  DBCD_Data <- read.csv(file_path)
  DBCD_Data <- na.omit(DBCD_Data)
  DBCD_Data$status <- as.integer(DBCD_Data$status)
  DBCD_Data$time <- as.numeric(DBCD_Data$time)

  X <- DBCD_Data[, -c(1, 2), drop = FALSE]
  X <- standardize_x(X)

  list(name = "DBCD", X = X, time = DBCD_Data$time, status = DBCD_Data$status)
}

load_chop_data <- function(file_path = "chop.csv") {
  if (!file.exists(file_path)) {
    stop("Cannot find ", file_path, ". Put chop.csv in the project working directory.")
  }

  Datachop <- read.csv(file_path)
  Datachop <- Datachop[, -1, drop = FALSE]
  Datachop <- na.omit(Datachop)
  Datachop$status <- as.integer(Datachop$status)
  Datachop$survtime <- as.numeric(Datachop$survtime)

  keep <- is.finite(Datachop$survtime) &
    Datachop$survtime > 0 &
    is.finite(Datachop$status)

  Datachop <- Datachop[keep, , drop = FALSE]

  X <- Datachop[, -c(1, 2), drop = FALSE]
  X <- standardize_x(X)

  list(name = "CHOP", X = X, time = Datachop$survtime, status = Datachop$status)
}

load_nki70 <- function() {
  data("nki70", package = "penalized")
  time <- as.numeric(nki70$time)
  status <- as.integer(nki70$event)
  X <- nki70[, !(names(nki70) %in% c("time", "event")), drop = FALSE]
  X <- standardize_x(X)

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
  X <- standardize_x(X)

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
  
  X <- standardize_x(X)
  
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
  load_pbc_data(),
  load_dbcd_data("DBCD_Data.csv"),
  load_chop_data("chop.csv"),
  load_nki70(),
  load_vdx(),
  load_ovarian()
)

real_results <- list()
auc_results <- list()
counter <- 1

for (dat in real_datasets) {
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

  if (ncol(dat$X) <= 100) {
    rsf_screen <- list(
      X_selected = dat$X,
      selected_variables = colnames(dat$X),
      importance = rep(NA_real_, ncol(dat$X)),
      threshold_rule = "No screening: p <= 100",
      threshold_value = NA_real_,
      original_p = ncol(dat$X),
      selected_p = ncol(dat$X)
    )
  } else {
    rsf_screen <- screen_by_rsf_quantile(
      X = dat$X,
      time = dat$time,
      status = dat$status,
      quantile_cutoff = rsf_quantile_cutoff,
      min_p = rsf_min_p,
      max_p = rsf_max_p,
      ntree = rsf_ntree,
      seed = 2026
    )
  }

  dat$X <- rsf_screen$X_selected

  cat("\nAfter RSF screening:\n")
  cat("Selected p =", ncol(dat$X), "\n")
  cat("Threshold rule =", rsf_screen$threshold_rule, "\n")
  cat("Threshold value =", rsf_screen$threshold_value, "\n")

  write.csv(
    data.frame(
      Variable = rsf_screen$selected_variables,
      Importance = as.numeric(rsf_screen$importance)
    ),
    paste0("results/real_data/RSF_selected_variables_", safe_dataset_name(dat$name), ".csv"),
    row.names = FALSE
  )

  split <- make_split(
    X = dat$X,
    time = dat$time,
    status = dat$status,
    train_fraction = train_fraction,
    seed = 2026
  )

  times_auc <- as.numeric(
    quantile(
      split$time_test[split$status_test == 1],
      probs = c(0.25, 0.50, 0.75),
      na.rm = TRUE
    )
  )
  times_auc <- unique(times_auc[is.finite(times_auc)])

  if (length(times_auc) < 2) {
    times_auc <- as.numeric(
      quantile(
        split$time_test,
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
    block_size = vi_block_size
  )
  time_agp <- proc.time()[3] - t0

  t0 <- proc.time()[3]
  fit_gl <- fit_gl_blended_igh_vi(
    split$X_train,
    split$time_train,
    split$status_train,
    max_iter = vi_max_iter,
    block_size = vi_block_size
  )
  time_gl <- proc.time()[3] - t0

  fit_lasso$method <- "Cox-Lasso"
  fit_enet$method <- "Cox-ElasticNet"
  fit_bayes$method <- "Bayesian-Lasso-Cox"
  fit_agp$method <- "AGP-VI"
  fit_gl$method <- "Global-Local Gaussian–IGH-VI"

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
    result$RSF_Selected_p <- ncol(dat$X)
    result$RSF_Threshold_Rule <- rsf_screen$threshold_rule
    result$RSF_Threshold_Value <- as.numeric(rsf_screen$threshold_value)
    result$Time_s <- as.numeric(unname(runtimes[method_name]))

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
        "Time_s"
      )
    ]

    real_results[[counter]] <- result
    counter <- counter + 1
  }

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
    rank_cindex = rank(-C_index_num, ties.method = "min"),
    rank_brier = rank(Brier_Score_num, ties.method = "min"),
    rank_auc = rank(-AUC_num, ties.method = "min"),
    overall_rank_score = rank_cindex + rank_brier + rank_auc,
    Best_Method = overall_rank_score == min(overall_rank_score, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    C_index = ifelse(Best_Method, paste0("\\textbf{", sprintf("%.4f", C_index_num), "}"), sprintf("%.4f", C_index_num)),
    Brier_Score = ifelse(Best_Method, paste0("\\textbf{", sprintf("%.4f", Brier_Score_num), "}"), sprintf("%.4f", Brier_Score_num)),
    Mean_Time_Dependent_AUC = ifelse(Best_Method, paste0("\\textbf{", sprintf("%.4f", AUC_num), "}"), sprintf("%.4f", AUC_num)),
    Mean_Selected = ifelse(Best_Method, paste0("\\textbf{", Mean_Selected_num, "}"), as.character(Mean_Selected_num))
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
  "Dataset & Method & C-index & Brier Score & Mean Time-Dependent AUC & Mean Selected \\\\",
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
