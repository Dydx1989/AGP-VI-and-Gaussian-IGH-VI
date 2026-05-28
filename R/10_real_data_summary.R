evaluate_real_model <- function(fit, X_train, time_train, status_train, X_test, time_test, status_test, threshold = 0.20, z_threshold = 2.58, selection_rule = "and") {
  result <- evaluate_model(fit, X_train, time_train, status_train, X_test, time_test, status_test, beta_true = NULL, threshold = threshold, z_threshold = z_threshold, selection_rule = selection_rule)
  selection <- select_variables(fit, threshold, z_threshold, selection_rule)
  result$Selected <- selection$n_selected
  result
}

summarise_real_results <- function(raw_results) {
  raw_results |>
    dplyr::group_by(Dataset, Method) |>
    dplyr::summarise(
      C_index = mean(C_index, na.rm = TRUE),
      SD_C_index = sd(C_index, na.rm = TRUE),
      Brier_Score = mean(Brier_Score, na.rm = TRUE),
      SD_Brier = sd(Brier_Score, na.rm = TRUE),
      Selected = mean(Selected, na.rm = TRUE),
      SD_Selected = sd(Selected, na.rm = TRUE),
      Time_s = mean(Time_s, na.rm = TRUE),
      SD_Time_s = sd(Time_s, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(Dataset, factor(Method, levels = c("Cox-Lasso", "Horseshoe-Cox", "AGP-VI", "Gaussian--IGH-VI")))
}

write_real_latex_table <- function(summary_results, file) {
  lines <- c("\\renewcommand{\\arraystretch}{1.5}", "\\begin{table}[ht]", "\\tiny", "\\caption{Real-Data Performance of Cox-Lasso, Horseshoe-Cox, AGP-VI, and Gaussian--IGH-VI}", "\\resizebox{\\textwidth}{!}{", "\\begin{tabular}{|l|l|c|c|c|c|}", "\\hline", "Dataset & Method & C-index & Brier Score & Selected Variables & Time (s) \\\", "\\hline")
  for (dataset_name in unique(summary_results$Dataset)) {
    sub <- summary_results[summary_results$Dataset == dataset_name, ]
    best_c <- max(sub$C_index, na.rm = TRUE)
    best_b <- min(sub$Brier_Score, na.rm = TRUE)
    best_t <- min(sub$Time_s, na.rm = TRUE)
    for (i in seq_len(nrow(sub))) {
      row <- sub[i, ]
      dataset_cell <- ifelse(i == 1, dataset_name, "")
      c_value <- fmt(row$C_index, 4)
      b_value <- fmt(row$Brier_Score, 6)
      selected_value <- fmt(row$Selected, 2)
      time_value <- fmt(row$Time_s, 2)
      if (!is.na(row$C_index) && abs(row$C_index - best_c) < 1e-12) c_value <- paste0("\\textbf{", c_value, "}")
      if (!is.na(row$Brier_Score) && abs(row$Brier_Score - best_b) < 1e-12) b_value <- paste0("\\textbf{", b_value, "}")
      if (!is.na(row$Time_s) && abs(row$Time_s - best_t) < 1e-12) time_value <- paste0("\\textbf{", time_value, "}")
      lines <- c(lines, paste0(paste(dataset_cell, as.character(row$Method), c_value, b_value, selected_value, time_value, sep = " & "), " \\\"))
    }
    lines <- c(lines, "\\hline")
  }
  lines <- c(lines, "\\end{tabular}", "}", "\\label{tab:real_data_performance}", "\\end{table}")
  writeLines(lines, con = file)
}

plot_real_data_comparison <- function(summary_results, save_path = NULL) {
  plot_data <- tidyr::pivot_longer(summary_results, cols = c("C_index", "Brier_Score", "Selected", "Time_s"), names_to = "Metric", values_to = "Value")
  plot_data$Metric <- dplyr::recode(plot_data$Metric, C_index = "C-index", Brier_Score = "Brier Score", Selected = "Selected Variables", Time_s = "Computation Time")
  plot_object <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dataset, y = Value, fill = Method)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~ Metric, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Real-Data Performance Across Datasets", x = "Dataset", y = "Value")
  if (!is.null(save_path)) ggplot2::ggsave(save_path, plot_object, width = 11, height = 7, dpi = 300)
  plot_object
}
