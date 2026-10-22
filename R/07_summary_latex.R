# ============================================================================
# Summary table and LaTeX export
# ============================================================================

summarise_results <- function(raw_results) {
  mean_safe <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    mean(x, na.rm = TRUE)
  }
  
  sd_safe <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    sd(x, na.rm = TRUE)
  }
  
  methods_order <- c(
    "Cox-Lasso",
    "Cox-ElasticNet",
    "Bayesian-Lasso-Cox",
    "AGP-VI",
    "Global-Local Gaussian–IGH-VI"
  )
  
  raw_results$TPR <- raw_results$TP / raw_results$n_signal
  
  raw_results$FPR <- raw_results$FP /
    pmax(raw_results$p - raw_results$n_signal, 1)
  
  output <- raw_results |>
    dplyr::group_by(p, Method) |>
    dplyr::summarise(
      C_index = mean_safe(C_index),
      Brier_Score = mean_safe(Brier_Score),
      Mean_TPR = mean_safe(TPR),
      SD_TPR = sd_safe(TPR),
      Mean_FPR = mean_safe(FPR),
      SD_FPR = sd_safe(FPR),
      Time_s = mean_safe(Time_s),
      .groups = "drop"
    )
  
  output$Method <- factor(
    output$Method,
    levels = methods_order
  )
  
  output <- output[
    order(output$p, output$Method),
  ]
  
  as.data.frame(output)
}


fmt <- function(x,
                digits = 4) {
  ifelse(
    is.na(x),
    "NA",
    formatC(
      x,
      format = "f",
      digits = digits
    )
  )
}


fmt_brier <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      abs(x) < 1e-4,
      paste0(
        "$",
        formatC(
          x,
          format = "e",
          digits = 2
        ),
        "$"
      ),
      formatC(
        x,
        format = "f",
        digits = 6
      )
    )
  )
}


bold_if_best <- function(value,
                         best_value,
                         formatted_value) {
  if (is.na(value) || is.na(best_value)) {
    return(formatted_value)
  }
  
  if (abs(value - best_value) < 1e-12) {
    return(
      paste0(
        "\\textbf{",
        formatted_value,
        "}"
      )
    )
  }
  
  formatted_value
}


write_latex_performance_table <- function(tab,
                                          file = "results/performance_table.tex") {
  if (!dir.exists(dirname(file))) {
    dir.create(
      dirname(file),
      recursive = TRUE
    )
  }
  
  lines <- c(
    "\\renewcommand{\\arraystretch}{1.7}",
    "\\begin{table}[ht]",
    "\\tiny",
    "\\caption{Performance Metrics and Computation Times for Cox-Lasso, Horseshoe-Cox, AGP-VI, and Gaussian--IGH-VI at Varying Predictor Dimensions}",
    "\\resizebox{\\textwidth}{!}{",
    "\\begin{tabular}{|c|l|c|c|c|c|c|c|c|}",
    "\\hline",
    "$p$ & Method & C-index & Brier Score & Mean TPR & SD TPR & Mean FPR & SD FPR & Time (s) \\\\",
    "\\hline"
  )
  
  for (pp in unique(tab$p)) {
    sub <- tab[
      tab$p == pp,
    ]
    
    best_c <- max(
      sub$C_index,
      na.rm = TRUE
    )
    
    best_brier <- min(
      sub$Brier_Score,
      na.rm = TRUE
    )
    
    best_tpr <- max(
      sub$Mean_TPR,
      na.rm = TRUE
    )
    
    best_fpr <- min(
      sub$Mean_FPR,
      na.rm = TRUE
    )
    
    best_time <- min(
      sub$Time_s,
      na.rm = TRUE
    )
    
    for (i in seq_len(nrow(sub))) {
      row <- sub[i, ]
      
      pcell <- ifelse(
        i == 1,
        as.character(pp),
        ""
      )
      
      c_value <- fmt(
        row$C_index,
        4
      )
      
      brier_value <- fmt_brier(
        row$Brier_Score
      )
      
      tpr_value <- fmt(
        row$Mean_TPR,
        4
      )
      
      sd_tpr_value <- fmt(
        row$SD_TPR,
        4
      )
      
      fpr_value <- fmt(
        row$Mean_FPR,
        4
      )
      
      sd_fpr_value <- fmt(
        row$SD_FPR,
        4
      )
      
      time_value <- fmt(
        row$Time_s,
        2
      )
      
      c_value <- bold_if_best(
        row$C_index,
        best_c,
        c_value
      )
      
      brier_value <- bold_if_best(
        row$Brier_Score,
        best_brier,
        brier_value
      )
      
      tpr_value <- bold_if_best(
        row$Mean_TPR,
        best_tpr,
        tpr_value
      )
      
      fpr_value <- bold_if_best(
        row$Mean_FPR,
        best_fpr,
        fpr_value
      )
      
      time_value <- bold_if_best(
        row$Time_s,
        best_time,
        time_value
      )
      
      latex_row <- paste(
        pcell,
        as.character(row$Method),
        c_value,
        brier_value,
        tpr_value,
        sd_tpr_value,
        fpr_value,
        sd_fpr_value,
        time_value,
        sep = " & "
      )
      
      lines <- c(
        lines,
        paste0(
          latex_row,
          " \\\\"
        )
      )
    }
    
    lines <- c(
      lines,
      "\\hline"
    )
  }
  
  lines <- c(
    lines,
    "\\end{tabular}",
    "}",
    "\\label{tab:performance_metrics}",
    "\\end{table}"
  )
  
  writeLines(
    lines,
    con = file
  )
}