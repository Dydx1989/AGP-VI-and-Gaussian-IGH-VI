# ============================================================================
# Prediction utilities
# ============================================================================

predict_risk <- function(X, beta) {
  
  beta[!is.finite(beta)] <- 0
  
  drop(as.matrix(X) %*% beta)
}

# ============================================================================
# Harrell C-index
# ============================================================================

harrell_cindex <- function(time, status, risk) {
  
  n <- length(time)
  
  comparable <- 0
  concordant <- 0
  tied <- 0
  
  for (i in 1:(n - 1)) {
    
    for (j in (i + 1):n) {
      
      if (time[i] == time[j]) {
        next
      }
      
      if (time[i] < time[j] && status[i] == 1) {
        
        comparable <- comparable + 1
        
        if (risk[i] > risk[j]) {
          concordant <- concordant + 1
        }
        
        if (risk[i] == risk[j]) {
          tied <- tied + 1
        }
      }
      
      if (time[j] < time[i] && status[j] == 1) {
        
        comparable <- comparable + 1
        
        if (risk[j] > risk[i]) {
          concordant <- concordant + 1
        }
        
        if (risk[j] == risk[i]) {
          tied <- tied + 1
        }
      }
    }
  }
  
  if (comparable == 0) {
    return(NA_real_)
  }
  
  (concordant + 0.5 * tied) / comparable
}

# ============================================================================
# Breslow baseline hazard
# ============================================================================

baseline_hazard_breslow <- function(X,
                                    time,
                                    status,
                                    beta) {
  
  eta <- drop(as.matrix(X) %*% beta)
  
  eta <- pmin(pmax(eta, -30), 30)
  
  event_times <- sort(unique(time[status == 1]))
  
  H0 <- numeric(length(event_times))
  
  cumhaz <- 0
  
  for (k in seq_along(event_times)) {
    
    tk <- event_times[k]
    
    dk <- sum(time == tk & status == 1)
    
    risk <- which(time >= tk)
    
    denominator <- sum(exp(eta[risk])) + 1e-12
    
    cumhaz <- cumhaz + dk / denominator
    
    H0[k] <- cumhaz
  }
  
  data.frame(
    time = event_times,
    H0 = H0
  )
}

# ============================================================================
# Survival prediction
# ============================================================================

predict_survival <- function(X,
                             beta,
                             basehaz,
                             eval_time) {
  
  eta <- drop(as.matrix(X) %*% beta)
  
  eta <- pmin(pmax(eta, -30), 30)
  
  H0t <- if (
    nrow(basehaz) == 0 ||
    all(basehaz$time > eval_time)
  ) {
    0
  } else {
    basehaz$H0[
      max(which(basehaz$time <= eval_time))
    ]
  }
  
  pmin(
    pmax(exp(-H0t * exp(eta)), 0),
    1
  )
}

# ============================================================================
# IPCW Brier score
# ============================================================================

ipcw_brier <- function(time,
                       status,
                       surv_prob,
                       eval_time) {
  
  censor_status <- 1 - status
  
  censor_fit <- survival::survfit(
    survival::Surv(time, censor_status) ~ 1
  )
  
  Ghat <- function(t) {
    
    pmax(
      summary(
        censor_fit,
        times = t,
        extend = TRUE
      )$surv,
      1e-6
    )
  }
  
  score <- numeric(length(time))
  
  for (i in seq_along(time)) {
    
    if (time[i] <= eval_time && status[i] == 1) {
      
      score[i] <- surv_prob[i]^2 / Ghat(time[i])
      
    } else if (time[i] > eval_time) {
      
      score[i] <- (1 - surv_prob[i])^2 / Ghat(eval_time)
      
    } else {
      
      score[i] <- 0
    }
  }
  
  mean(score)
}

# ============================================================================
# Variable selection
# ============================================================================

select_variables <- function(fit,
                             absolute_threshold = 0.15,
                             z_threshold = 2.33,
                             rule = c("or", "and"),
                             selection_method = c("auto", "credible_interval", "legacy", "nonzero"),
                             credible_level = 0.95,
                             nonzero_tol = 1e-8) {
  rule <- match.arg(rule)
  selection_method <- match.arg(selection_method)
  beta <- as.numeric(fit$beta)
  beta[!is.finite(beta)] <- 0
  p <- length(beta)

  has_cov <- !is.null(fit$covariance) &&
    is.matrix(fit$covariance) &&
    all(dim(fit$covariance) == c(p, p))

  method_used <- selection_method
  if (selection_method == "auto") {
    method_used <- if (has_cov) "credible_interval" else "nonzero"
  }

  se <- rep(NA_real_, p)
  z_score <- rep(NA_real_, p)
  ci_lower <- rep(NA_real_, p)
  ci_upper <- rep(NA_real_, p)

  if (method_used == "credible_interval") {
    if (!has_cov) {
      warning("credible_interval selection requested but covariance is unavailable; using nonzero rule.")
      method_used <- "nonzero"
    } else {
      se <- sqrt(pmax(diag(fit$covariance), 1e-12))
      z_score <- abs(beta) / pmax(se, 1e-12)
      zcrit <- stats::qnorm(1 - (1 - credible_level) / 2)
      ci_lower <- beta - zcrit * se
      ci_upper <- beta + zcrit * se
      selected <- (ci_lower > 0) | (ci_upper < 0)
    }
  }

  if (method_used == "legacy") {
    if (has_cov) {
      se <- sqrt(pmax(diag(fit$covariance), 1e-12))
      z_score <- abs(beta) / pmax(se, 1e-12)
      if (rule == "and") {
        selected <- abs(beta) > absolute_threshold & z_score > z_threshold
      } else {
        selected <- abs(beta) > absolute_threshold | z_score > z_threshold
      }
    } else {
      selected <- abs(beta) > absolute_threshold
    }
  }

  if (method_used == "nonzero") {
    selected <- abs(beta) > nonzero_tol
  }

  list(
    selected = as.logical(selected),
    beta = beta,
    se = se,
    z_score = z_score,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    n_selected = sum(selected, na.rm = TRUE),
    selection_method = method_used,
    credible_level = if (method_used == "credible_interval") credible_level else NA_real_
  )
}

# ============================================================================
# Selection metrics
# ============================================================================

selection_metrics <- function(fit,
                              beta_true,
                              threshold = 0.15,
                              z_threshold = 2.33,
                              rule = "or",
                              selection_method = "auto",
                              credible_level = 0.95) {
  
  selection <- select_variables(
    fit = fit,
    absolute_threshold = threshold,
    z_threshold = z_threshold,
    rule = rule,
    selection_method = selection_method,
    credible_level = credible_level
  )
  
  signal <- beta_true != 0
  
  c(
    TP = sum(selection$selected & signal),
    FP = sum(selection$selected & !signal)
  )
}

# ============================================================================
# Main evaluation
# ============================================================================

evaluate_model <- function(fit,
                           X_train,
                           time_train,
                           status_train,
                           X_test,
                           time_test,
                           status_test,
                           beta_true = NULL,
                           eval_time = NULL,
                           threshold = 0.15,
                           z_threshold = 2.33,
                           selection_rule = "or",
                           selection_method = "auto",
                           credible_level = 0.95) {
  
  beta_hat <- fit$beta
  
  if (
    length(beta_hat) != ncol(X_test) ||
    all(is.na(beta_hat))
  ) {
    
    return(
      data.frame(
        Method = fit$method,
        C_index = NA_real_,
        Brier_Score = NA_real_,
        TP = NA_real_,
        FP = NA_real_
      )
    )
  }
  
  beta_hat[!is.finite(beta_hat)] <- 0
  
  risk <- predict_risk(X_test, beta_hat)
  
  if (is.null(eval_time) || !is.finite(eval_time)) {
    
    eval_time <- if (
      sum(status_train == 1) > 0
    ) {
      
      as.numeric(
        stats::quantile(
          time_train[status_train == 1],
          0.50,
          na.rm = TRUE
        )
      )
      
    } else {
      
      as.numeric(
        stats::median(
          time_train,
          na.rm = TRUE
        )
      )
    }
  }
  
  basehaz <- baseline_hazard_breslow(
    X_train,
    time_train,
    status_train,
    beta_hat
  )
  
  surv <- predict_survival(
    X_test,
    beta_hat,
    basehaz,
    eval_time
  )
  
  if (!is.null(beta_true)) {
    
    sel <- selection_metrics(
      fit = fit,
      beta_true = beta_true,
      threshold = threshold,
      z_threshold = z_threshold,
      rule = selection_rule,
      selection_method = selection_method,
      credible_level = credible_level
    )
    
    TP <- sel["TP"]
    
    FP <- sel["FP"]
    
  } else {
    
    TP <- NA_real_
    
    FP <- NA_real_
  }
  
  data.frame(
    Method = fit$method,
    C_index = harrell_cindex(
      time_test,
      status_test,
      risk
    ),
    Brier_Score = ipcw_brier(
      time_test,
      status_test,
      surv,
      eval_time
    ),
    TP = as.numeric(TP),
    FP = as.numeric(FP)
  )
}