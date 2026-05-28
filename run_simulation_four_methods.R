# ============================================================================
# 09_main_simulation.R
# Main simulation study
# ============================================================================

rm(list = ls())

# ----------------------------------------------------------------------------
# Load packages and source functions
# ----------------------------------------------------------------------------

source("R/00_packages.R")
load_required_packages()

source("R/01_simulate_data.R")
source("R/02_common_functions.R")
source("R/03_fit_agp_vi.R")
source("R/04b_fit_gl_blended_igh_vi.R")
source("R/05_fit_existing_methods.R")
source("R/06_evaluation.R")
source("R/07_summary_latex.R")
source("R/08_plots.R")

# ----------------------------------------------------------------------------
# Create folders
# ----------------------------------------------------------------------------

if (!dir.exists("results/simulation")) {
  dir.create("results/simulation", recursive = TRUE)
}

if (!dir.exists("figures/simulation")) {
  dir.create("figures/simulation", recursive = TRUE)
}

# ----------------------------------------------------------------------------
# Simulation settings
# ----------------------------------------------------------------------------

set.seed(202236)

n_rep <- 5
n <- 100

p_grid <- c(
  50,
  100,
  200,300,400,500
)

n_signal <- 10
rho <- 0.15
beta_value <- 1.20
censoring_rate <- 0.50
train_fraction <- 0.70
absolute_threshold <- 0.15
z_threshold <- 2.33
raw_list <- list()
counter <- 1
saved_vi_plots <- FALSE

# ----------------------------------------------------------------------------
# Main simulation loop
# ----------------------------------------------------------------------------

for (pp in p_grid) {
  
  cat("\n========================================\n")
  cat("Simulation for p =", pp, "\n")
  cat("========================================\n")
  
  vi_max_iter <- ifelse(pp >= 300, 250, 300)
  
  for (rr in seq_len(n_rep)) {
    
    cat("Replication", rr, "of", n_rep, "\n")
    
    # ------------------------------------------------------------------------
    # Generate data
    # ------------------------------------------------------------------------
    
    simulated_data <- simulate_cox_data(
      n = n,
      p = pp,
      n_signal = n_signal,
      rho = rho,
      beta_value = beta_value,
      censoring_rate = censoring_rate,
      seed = 100000 + pp * 100 + rr
    )
    
    split_data <- train_test_split(
      X = simulated_data$X,
      time = simulated_data$time,
      status = simulated_data$status,
      train_fraction = train_fraction,
      seed = 200000 + pp * 100 + rr
    )
    
    eval_time <- if (
      sum(split_data$status_train == 1) > 0
    ) {
      as.numeric(
        quantile(
          split_data$time_train[
            split_data$status_train == 1
          ],
          0.50
        )
      )
    } else {
      median(split_data$time_train)
    }
    
    # ------------------------------------------------------------------------
    # Cox-Lasso
    # ------------------------------------------------------------------------
    
    start_time <- proc.time()[3]
    
    fit_cox_lasso_object <- tryCatch(
      fit_cox_lasso(
        X = split_data$X_train,
        time = split_data$time_train,
        status = split_data$status_train
      ),
      error = function(e) {
        warning(
          "Cox-Lasso failed: ",
          conditionMessage(e)
        )
        
        list(
          method = "Cox-Lasso",
          beta = rep(NA_real_, pp)
        )
      }
    )
    
    time_cox_lasso <- proc.time()[3] - start_time
    
    result_cox_lasso <- evaluate_model(
      fit = fit_cox_lasso_object,
      X_train = split_data$X_train,
      time_train = split_data$time_train,
      status_train = split_data$status_train,
      X_test = split_data$X_test,
      time_test = split_data$time_test,
      status_test = split_data$status_test,
      beta_true = simulated_data$beta_true,
      eval_time = eval_time,
      threshold = absolute_threshold,
      z_threshold = z_threshold
    )
    
    result_cox_lasso$Time_s <- time_cox_lasso
    
    # ------------------------------------------------------------------------
    # Cox-ElasticNet
    # ------------------------------------------------------------------------
    
    start_time <- proc.time()[3]
    
    fit_elasticnet_object <- tryCatch(
      fit_cox_elasticnet(
        X = split_data$X_train,
        time = split_data$time_train,
        status = split_data$status_train,
        alpha = 0.5
      ),
      error = function(e) {
        warning("Cox-ElasticNet failed: ", conditionMessage(e))
        
        list(
          method = "Cox-ElasticNet",
          beta = rep(NA_real_, pp)
        )
      }
    )
    
    time_elasticnet <- proc.time()[3] - start_time
    
    result_elasticnet <- evaluate_model(
      fit = fit_elasticnet_object,
      X_train = split_data$X_train,
      time_train = split_data$time_train,
      status_train = split_data$status_train,
      X_test = split_data$X_test,
      time_test = split_data$time_test,
      status_test = split_data$status_test,
      beta_true = simulated_data$beta_true,
      eval_time = eval_time,
      threshold = absolute_threshold,
      z_threshold = z_threshold
    )
    
    result_elasticnet$Method <- "Cox-ElasticNet"
    result_elasticnet$Time_s <- time_elasticnet
    
    # ------------------------------------------------------------------------
    # Bayesian Lasso Cox
    # ------------------------------------------------------------------------
    
    start_time <- proc.time()[3]
    
    fit_bayes_lasso_object <- tryCatch(
      fit_bayesian_lasso_cox(
        X = split_data$X_train,
        time = split_data$time_train,
        status = split_data$status_train,
        chains = 2,
        iter = 500,
        warmup = 250,
        seed = 400000 + pp * 100 + rr
      ),
      error = function(e) {
        warning("Bayesian Lasso Cox failed: ", conditionMessage(e))
        
        list(
          method = "Bayesian-Lasso-Cox",
          beta = rep(NA_real_, pp)
        )
      }
    )
    
    time_bayes_lasso <- proc.time()[3] - start_time
    
    result_bayes_lasso <- evaluate_model(
      fit = fit_bayes_lasso_object,
      X_train = split_data$X_train,
      time_train = split_data$time_train,
      status_train = split_data$status_train,
      X_test = split_data$X_test,
      time_test = split_data$time_test,
      status_test = split_data$status_test,
      beta_true = simulated_data$beta_true,
      eval_time = eval_time,
      threshold = absolute_threshold,
      z_threshold = z_threshold
    )
    
    result_bayes_lasso$Method <- "Bayesian-Lasso-Cox"
    result_bayes_lasso$Time_s <- time_bayes_lasso
    # ------------------------------------------------------------------------
    # AGP-VI
    # ------------------------------------------------------------------------
    
    start_time <- proc.time()[3]
    
    fit_agp_object <- tryCatch(
      fit_agp_vi(
        X = split_data$X_train,
        time = split_data$time_train,
        status = split_data$status_train,
        max_iter = vi_max_iter,
        eta = 0.001,
        tol = 1e-4,
        init_var = 0.5,
        min_var = 1e-4,
        max_var = 5,
        ridge = 0.05,
        grad_clip = 5,
        block_size = 5,
        corr_shrinkage = 0.3,
        verbose = FALSE
      ),
      error = function(e) {
        
        warning(
          "AGP-VI failed: ",
          conditionMessage(e)
        )
        
        list(
          method = "AGP-VI",
          beta = rep(NA_real_, pp),
          covariance = diag(NA_real_, pp),
          elbo = numeric()
        )
      }
    )
    
    time_agp <- proc.time()[3] - start_time
    
    result_agp <- evaluate_model(
      fit = fit_agp_object,
      X_train = split_data$X_train,
      time_train = split_data$time_train,
      status_train = split_data$status_train,
      X_test = split_data$X_test,
      time_test = split_data$time_test,
      status_test = split_data$status_test,
      beta_true = simulated_data$beta_true,
      eval_time = eval_time,
      threshold = absolute_threshold,
      z_threshold = z_threshold
    )
    
    result_agp$Time_s <- time_agp
    
    
    # ------------------------------------------------------------------------
    # Global-Local Gaussian–IGH-VI
    # ------------------------------------------------------------------------
    
    start_time <- proc.time()[3]
    
    fit_gl_igh_object <- tryCatch(
      fit_gl_blended_igh_vi(
        X = split_data$X_train,
        time = split_data$time_train,
        status = split_data$status_train,
        a = 1,
        b = 0.5,
        a0 = 1,
        b0 = 2,
        blend_weight = 0.90,
        max_iter = vi_max_iter,
        eta = 0.001,
        tol = 1e-4,
        init_var = 0.5,
        min_var = 1e-4,
        max_var = 5,
        min_global_var = 0.05,
        max_global_var = 10,
        ridge = 0.04,
        grad_clip = 5,
        block_size = 10,
        corr_shrinkage = 0.3,
        verbose = FALSE
      ),
      error = function(e) {
        warning(
          "Global-Local Gaussian–IGH-VI failed: ",
          conditionMessage(e)
        )
        
        list(
          method = "Global-Local Gaussian–IGH-VI",
          beta = rep(NA_real_, pp),
          covariance = diag(NA_real_, pp),
          tau2 = rep(NA_real_, pp),
          tau0_2 = NA_real_,
          effective_var = rep(NA_real_, pp),
          elbo = numeric()
        )
      }
    )
    
    if (is.null(fit_gl_igh_object$method) || is.na(fit_gl_igh_object$method)) {
      fit_gl_igh_object$method <- "Global-Local Gaussian–IGH-VI"
    }
    
    time_gl_igh <- proc.time()[3] - start_time
    
    fit_gl_igh_object$method <- "Global-Local Gaussian–IGH-VI"
    
    result_gl_igh <- evaluate_model(
      fit = fit_gl_igh_object,
      X_train = split_data$X_train,
      time_train = split_data$time_train,
      status_train = split_data$status_train,
      X_test = split_data$X_test,
      time_test = split_data$time_test,
      status_test = split_data$status_test,
      beta_true = simulated_data$beta_true,
      eval_time = eval_time,
      threshold = absolute_threshold,
      z_threshold = z_threshold
    )
    
    result_gl_igh$Method <- "Global-Local Gaussian–IGH-VI"
    result_gl_igh$Time_s <- time_gl_igh
    
    # ------------------------------------------------------------------------
    # Save plots once
    # ------------------------------------------------------------------------
    
    if (!saved_vi_plots) {
      
      # --------------------------------------------------------------------
      # ELBO plots
      # --------------------------------------------------------------------
      
      plot_vi_elbo(
        fit_agp_object,
        "AGP-VI",
        "figures/simulation/AGP_VI_ELBO_convergence.png"
      )
      
      plot_vi_elbo(
        fit_igh_object,
        "Global-Local Gaussian–IGH-VI",
        "figures/simulation/Gaussian_IGH_VI_ELBO_convergence.png"
      )
      
      plot_vi_elbo(
        fit_gl_igh_object,
        "Improved Tempered GL-Blended Gaussian--IGH-VI",
        "figures/simulation/GL_IGH_VI_ELBO_convergence.png"
      )
      
      # --------------------------------------------------------------------
      # Coefficient recovery across all methods
      # --------------------------------------------------------------------
      
      fit_list <- list(
        "Cox-Lasso" = fit_cox_lasso_object,
        "Cox-ElasticNet" = fit_elasticnet_object,
        "Bayesian-Lasso-Cox" = fit_bayes_lasso_object,
        "AGP-VI" = fit_agp_object,
        "Global-Local Gaussian–IGH-VI" = fit_gl_igh_object
      )
      
      plot_coefficient_recovery_all(
        beta_true = simulated_data$beta_true,
        fit_list = fit_list,
        save_path = "figures/simulation/coefficient_recovery_all_methods.png"
      )
      
      saved_vi_plots <- TRUE
    }
    # ------------------------------------------------------------------------
    # Combine results
    # ------------------------------------------------------------------------
    
    combined_result <- rbind(
      result_cox_lasso,
      result_elasticnet,
      result_bayes_lasso,
      result_agp,
      result_gl_igh
    )
    
    
    combined_result$p <- pp
    combined_result$n <- n
    combined_result$n_signal <- n_signal
    combined_result$rho <- rho
    combined_result$censoring_target <- censoring_rate
    combined_result$censoring_observed <- simulated_data$observed_censoring
    combined_result$replication <- rr
    
    raw_list[[counter]] <- combined_result
    
    counter <- counter + 1
    
    partial_results <- do.call(
      rbind,
      raw_list
    )
    
    write.csv(
      partial_results,
      "results/simulation/raw_results_partial.csv",
      row.names = FALSE
    )
  }
}

# ----------------------------------------------------------------------------
# Final results
# ----------------------------------------------------------------------------

raw_results <- do.call(
  rbind,
  raw_list
)

write.csv(
  raw_results,
  "results/simulation/raw_results.csv",
  row.names = FALSE
)

  summary_results <- summarise_results(
  raw_results
)

write.csv(
  summary_results,
  "results/simulation/summary_results.csv",
  row.names = FALSE
)

write_latex_performance_table(
  summary_results,
  "results/simulation/performance_table.tex"
)

# ----------------------------------------------------------------------------
# Final plots
# ----------------------------------------------------------------------------

plot_tpr_fpr_all(
  summary_results,
  "figures/simulation/TPR_FPR_all_methods.png"
)

plot_prediction_performance(
  summary_results,
  "figures/simulation/prediction_performance_all_methods.png"
)

plot_integrated_brier_score(
  summary_results,
  "figures/simulation/integrated_brier_score.png"
)

plot_runtime_scaling(
  summary_results,
  "figures/simulation/runtime_scaling_all_methods.png"
)

plot_simulation_boxplots(
  raw_results,
  "figures/simulation/simulation_boxplots_all_methods.png"
)
# ----------------------------------------------------------------------------
# Session info
# ----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = "results/simulation/session_info.txt"
)

cat(
  "\nSimulation complete. Files saved in results/simulation/ and figures/simulation/.\n"
)

print(summary_results)

