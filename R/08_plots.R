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
  if (is.null(fit$elbo) || length(fit$elbo) < 2) {
    warning("No usable ELBO trace for ", method_name)
    return(NULL)
  }
  
  elbo <- as.numeric(fit$elbo)
  
  plot_data <- data.frame(
    Iteration = seq_along(elbo),
    ELBO = elbo
  )
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(Iteration, ELBO)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste(method_name, "ELBO Convergence"),
      x = "Iteration",
      y = "ELBO"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  
  if (!is.null(save_path)) {
    
    dir.create(
      dirname(save_path),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    # ========================================================
    # DIAGNOSTIC: exact ELBO object being saved
    # ========================================================
    
    cat("\n===== EXACT OBJECT BEING SAVED =====\n")
    cat("Method:", method_name, "\n")
    cat(
      "Status:",
      ifelse(
        is.null(fit$convergence_status),
        "unknown",
        fit$convergence_status
      ),
      "\n"
    )
    cat("Iterations:", length(elbo), "\n")
    cat("First ELBO:", head(elbo, 1), "\n")
    cat("Final ELBO:", tail(elbo, 1), "\n")
    cat("Save path:", normalizePath(save_path, mustWork = FALSE), "\n")
    
    cat("\nFirst 5 ELBO values:\n")
    print(head(elbo, 5))
    
    cat("\nLast 5 ELBO values:\n")
    print(tail(elbo, 5))
    
    # Remove old figure before saving
    if (file.exists(save_path)) {
      file.remove(save_path)
    }
    
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 7,
      height = 5,
      dpi = 300
    )
    
    cat(
      "[FIGURE SAVED]",
      method_name,
      "| iterations =", length(elbo),
      "| final ELBO =", tail(elbo, 1),
      "|", normalizePath(save_path, mustWork = FALSE),
      "\n"
    )
  }
  
  p
}

# ============================================================================
# Rigorous ELBO convergence diagnostics: tail relative change against the
# algorithm's own stopping tolerance, mirroring exactly what the fitting
# functions' internal convergence check computes
# (rel_change <- abs(elbo[t] - elbo[t-1]) / (abs(elbo[t-1]) + 1e-8)). A small
# max(|rel_change|) over the tail is stronger evidence than a plot "looking
# flat": it is literally the same criterion the algorithm itself uses to
# decide convergence, checked independently and over a whole window rather
# than trusting a single accepted stopping point. This does NOT by itself
# rule out a frozen (all-fallback) tail, since a no-op fallback also gives
# rel_change = 0 trivially -- always report n_fallback_used alongside this,
# not instead of it (see PATCH_NOTES_GL_IGH_VI_oscillation.md, Attempt 6).
# ============================================================================
elbo_convergence_diagnostics <- function(elbo_trace, n_tail = 20, tol = 1e-4) {
  n <- length(elbo_trace)
  if (n < 3) {
    return(list(
      n_iterations = n,
      tail_elbo = elbo_trace,
      diff_tail = NA_real_,
      max_abs_diff_tail = NA_real_,
      rel_change_tail = NA_real_,
      max_abs_tail_rel_change = NA_real_,
      all_below_tol = NA
    ))
  }
  n_tail_actual <- min(n_tail, n - 1)
  tail_elbo <- utils::tail(elbo_trace, n_tail_actual + 1)
  diff_tail <- diff(tail_elbo)
  rel_change_full <- abs(diff(elbo_trace)) / (abs(head(elbo_trace, -1)) + 1e-8)
  rel_change_tail <- utils::tail(rel_change_full, n_tail_actual)
  list(
    n_iterations = n,
    tail_elbo = tail_elbo,
    diff_tail = diff_tail,
    max_abs_diff_tail = max(abs(diff_tail)),
    rel_change_tail = rel_change_tail,
    max_abs_tail_rel_change = max(rel_change_tail),
    all_below_tol = all(rel_change_tail < tol)
  )
}

# Prints a one-line summary plus the full diagnostic vectors (so they can be
# inspected directly, matching the tail(elbo_trace, 20) / diff / rel_change
# workflow), and returns the diagnostics list invisibly for programmatic use.
print_elbo_diagnostics <- function(fit, method_name, dataset_name, tol = 1e-4, n_tail = 20) {
  if (is.null(fit$elbo) || length(fit$elbo) < 3) {
    cat(sprintf(
      "  [ELBO diagnostic] %s %s: too few iterations for tail diagnostics (n=%d)\n",
      dataset_name, method_name, length(fit$elbo)
    ))
    return(invisible(NULL))
  }
  # The optimizer declares convergence after `patience` consecutive accepted
  # updates below tolerance.  Checking all of the last 20 changes can therefore
  # misleadingly print FALSE even for a legitimately converged fit (because
  # earlier values in that 20-step window may predate the convergence streak).
  # Report both the broad tail diagnostic and the actual stopping-window check.
  d <- elbo_convergence_diagnostics(fit$elbo, n_tail = n_tail, tol = tol)
  n_fb <- if (!is.null(fit$n_fallback_used)) fit$n_fallback_used else NA_integer_
  patience_used <- if (!is.null(fit$patience)) fit$patience else 5L
  stop_window <- utils::tail(d$rel_change_tail, min(patience_used, length(d$rel_change_tail)))
  stop_window_all_below <- length(stop_window) > 0 && all(stop_window < tol)
  cat(sprintf(
    "  [ELBO diagnostic] %s %s: n_iter=%d, n_fallback_used=%s, max|rel_change| over last %d = %.3e; stopping-window(last %d) all below tol(%.0e) = %s\n",
    dataset_name, method_name, d$n_iterations, as.character(n_fb), n_tail,
    d$max_abs_tail_rel_change, length(stop_window), tol, stop_window_all_below
  ))
  cat("    tail(elbo_trace, ", n_tail + 1, "):\n", sep = "")
  print(round(d$tail_elbo, 4))
  cat("    diff(tail):\n")
  print(round(d$diff_tail, 6))
  cat("    tail(rel_change, ", n_tail, "):\n", sep = "")
  print(signif(d$rel_change_tail, 4))
  invisible(d)
}


# ============================================================================
# Two-panel ELBO convergence plot for publication. Panel A shows the full
# trace (so the initial rapid climb from a poor starting point remains
# visible). Panel B shows only the tail (default: the last `zoom_last_n`
# iterations, or the whole trace if shorter), plotted at its own y-axis
# scale, so genuine late-stage stabilization -- or the lack of it -- is
# visible rather than being compressed flat by the much larger scale of the
# initial climb. Requires gridExtra (already loaded by R/00_packages.R).
# ============================================================================
plot_vi_elbo_two_panel <- function(fit, method_name, save_path = NULL, zoom_last_n = 350) {
  if (is.null(fit$elbo) || length(fit$elbo) < 2) return(NULL)
  
  n <- length(fit$elbo)
  plot_data <- data.frame(Iteration = seq_len(n), ELBO = fit$elbo)
  zoom_start <- max(1, n - zoom_last_n + 1)
  zoom_data <- plot_data[plot_data$Iteration >= zoom_start, ]
  
  panel_a <- ggplot2::ggplot(plot_data, ggplot2::aes(Iteration, ELBO)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.0) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste(method_name, "-- full trace (Panel A)"),
      x = "Iteration", y = "ELBO"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10))
  
  panel_b <- ggplot2::ggplot(zoom_data, ggplot2::aes(Iteration, ELBO)) +
    ggplot2::geom_line(linewidth = 1, color = "darkred") +
    ggplot2::geom_point(size = 1.0, color = "darkred") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0(method_name, " -- zoomed tail, iterations ", zoom_start, "-", n, " (Panel B)"),
      x = "Iteration", y = "ELBO"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10))
  
  combined <- gridExtra::grid.arrange(panel_a, panel_b, ncol = 1)
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, combined, width = 7, height = 8, dpi = 300)
  }
  
  invisible(combined)
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