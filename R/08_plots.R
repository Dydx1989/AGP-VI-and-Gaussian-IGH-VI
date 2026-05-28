# ============================================================================
# 08_plots.R
# Plot functions for all methods + time-dependent Brier score
# ============================================================================

method_levels <- c(
  "Cox-Lasso",
  "Cox-ElasticNet",
  "Bayesian-Lasso-Cox",
  "AGP-VI",
  "Global-Local Gaussian–IGH-VI"
)

plot_vi_elbo <- function(fit, method_name, save_path = NULL) {
  if (is.null(fit$elbo) || length(fit$elbo) < 2) return(NULL)
  
  plot_data <- data.frame(
    Iteration = seq_along(fit$elbo),
    ELBO = fit$elbo
  )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(Iteration, ELBO)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste(method_name, "ELBO Convergence"),
      x = "Iteration",
      y = "ELBO"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 7, height = 5, dpi = 300)
  }
  
  p
}


plot_coefficient_recovery_all <- function(beta_true,
                                          fit_list,
                                          save_path = NULL) {
  
  make_df <- function(fit, method_name) {
    if (is.null(fit$covariance)) {
      posterior_sd <- rep(NA_real_, length(fit$beta))
    } else {
      posterior_sd <- sqrt(pmax(diag(fit$covariance), 1e-12))
    }
    
    data.frame(
      Predictor = seq_along(beta_true),
      True_Beta = beta_true,
      Estimated_Beta = fit$beta,
      Lower = fit$beta - 1.96 * posterior_sd,
      Upper = fit$beta + 1.96 * posterior_sd,
      Method = method_name
    )
  }
  
  plot_data <- do.call(
    rbind,
    Map(make_df, fit_list, names(fit_list))
  )
  
  plot_data$Method <- factor(plot_data$Method, levels = method_levels)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Predictor)) +
    ggplot2::geom_point(
      ggplot2::aes(y = True_Beta),
      color = "red",
      size = 2,
      shape = 17
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = Estimated_Beta),
      size = 1.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Lower, ymax = Upper),
      width = 0.2,
      alpha = 0.40,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~ Method, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Coefficient Recovery Across Methods",
      subtitle = "Red triangles: true coefficients; points: estimated coefficients",
      x = "Predictor Index",
      y = "Coefficient Value"
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 12, height = 7, dpi = 300)
  }
  
  p
}


plot_tpr_fpr_all <- function(summary_results, save_path = NULL) {
  plot_data <- tidyr::pivot_longer(
    summary_results,
    cols = c("Mean_TPR", "Mean_FPR"),
    names_to = "Metric",
    values_to = "Value"
  )
  
  plot_data$Metric <- dplyr::recode(
    plot_data$Metric,
    Mean_TPR = "True Positive Rate",
    Mean_FPR = "False Positive Rate"
  )
  
  plot_data$Method <- factor(plot_data$Method, levels = method_levels)
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = factor(p), y = Value, fill = Method)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~ Metric) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Variable Selection Performance Across Methods",
      x = "Number of Predictors",
      y = "Proportion"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 11, height = 6, dpi = 300)
  }
  
  p
}


plot_prediction_performance <- function(summary_results, save_path = NULL) {
  plot_data <- tidyr::pivot_longer(
    summary_results,
    cols = c("C_index", "Brier_Score"),
    names_to = "Metric",
    values_to = "Value"
  )
  
  plot_data$Metric <- dplyr::recode(
    plot_data$Metric,
    C_index = "C-index",
    Brier_Score = "Brier Score"
  )
  
  plot_data$Method <- factor(plot_data$Method, levels = method_levels)
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = p, y = Value, color = Method, group = Method)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ Metric, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Prediction Performance Across Predictor Dimensions",
      x = "Number of Predictors",
      y = "Value"
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 10, height = 6, dpi = 300)
  }
  
  p
}


plot_runtime_scaling <- function(summary_results, save_path = NULL) {
  summary_results$Method <- factor(summary_results$Method, levels = method_levels)
  
  p <- ggplot2::ggplot(
    summary_results,
    ggplot2::aes(x = p, y = Time_s, color = Method, group = Method)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Computational Time Scaling",
      x = "Number of Predictors",
      y = "Mean Time (seconds, log scale)"
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 9, height = 5, dpi = 300)
  }
  
  p
}


plot_simulation_boxplots <- function(raw_results, save_path = NULL) {
  raw_results$TPR <- raw_results$TP / raw_results$n_signal
  raw_results$FPR <- raw_results$FP / (raw_results$p - raw_results$n_signal)
  
  long_data <- tidyr::pivot_longer(
    raw_results,
    cols = c("C_index", "Brier_Score", "TPR", "FPR", "Time_s"),
    names_to = "Metric",
    values_to = "Value"
  )
  
  long_data$Metric <- dplyr::recode(
    long_data$Metric,
    C_index = "C-index",
    Brier_Score = "Brier Score",
    TPR = "True Positive Rate",
    FPR = "False Positive Rate",
    Time_s = "Computation Time"
  )
  
  long_data$Method <- factor(long_data$Method, levels = method_levels)
  
  p <- ggplot2::ggplot(
    long_data,
    ggplot2::aes(x = Method, y = Value, fill = Method)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.4) +
    ggplot2::facet_grid(Metric ~ p, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Simulation Performance Across Methods and Predictor Dimensions",
      x = "",
      y = "Value"
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 15, height = 10, dpi = 300)
  }
  
  p
}


# ============================================================================
# Time-dependent Brier score and IBS plots
# ============================================================================

plot_time_dependent_brier <- function(brier_results, save_path = NULL) {
  brier_results$Method <- factor(brier_results$Method, levels = method_levels)
  
  p <- ggplot2::ggplot(
    brier_results,
    ggplot2::aes(
      x = Time,
      y = Brier_Score,
      color = Method,
      group = Method
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ p, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Time-dependent Brier Score",
      x = "Evaluation Time",
      y = "Brier Score"
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 11, height = 6, dpi = 300)
  }
  
  p
}


plot_integrated_brier_score <- function(ibs_results, save_path = NULL) {
  ibs_results$Method <- factor(ibs_results$Method, levels = method_levels)
  
  p <- ggplot2::ggplot(
    ibs_results,
    ggplot2::aes(
      x = p,
      y = Brier_Score,
      color = Method,
      group = Method
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Integrated Brier Score Across Predictor Dimensions",
      x = "Number of Predictors",
      y = "Integrated Brier Score"
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 9, height = 5, dpi = 300)
  }
  
  p
}