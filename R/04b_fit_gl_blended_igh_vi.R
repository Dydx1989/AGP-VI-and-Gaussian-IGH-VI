# ============================================================================
# Global-Local Gaussian--IGH-VI
# V10: curvature-aware variational mean update, coherent full-state
# backtracking, separate stalled/converged states, and no componentwise
# gradient clipping (which created artificial coefficient plateaus).
# ============================================================================

fit_gl_blended_igh_vi <- function(X,
                                  time,
                                  status,
                                  a = 1,
                                  b = 0.5,
                                  a0 = 1,
                                  b0 = 2,
                                  blend_weight = 0.90,
                                  max_iter = 300,
                                  eta = 1.0,
                                  tol = 1e-4,
                                  init_var = 0.5,
                                  min_var = 1e-4,
                                  max_var = 5,
                                  min_global_var = 0.05,
                                  max_global_var = 10,
                                  ridge = 0.04,
                                  block_size = 10,
                                  corr_shrinkage = 0.3,
                                  variance_damping = 0.05,
                                  max_direction_norm = 5,
                                  patience = 5,
                                  stall_patience = 20,
                                  max_backtrack = 20,
                                  rescue_mean_only = TRUE,
                                  rescue_raw_gradient = TRUE,
                                  rescue_max_backtrack = 40,
                                  min_rescue_gain_rel = 1e-12,
                                  m_init = NULL,
                                  S_init = NULL,
                                  tau2_init = NULL,
                                  tau0_2_init = NULL,
                                  verbose = FALSE) {
  X <- as.matrix(X)
  p <- ncol(X)
  if (p < 1) stop("X must contain at least one predictor.")

  blend_weight <- min(max(blend_weight, 0), 1)
  variance_damping <- min(max(variance_damping, 1e-3), 1)
  eta <- max(as.numeric(eta), 1e-8)

  if (!is.null(m_init)) {
    m <- as.numeric(m_init)
    if (length(m) != p) stop("m_init length does not match ncol(X)")
    m[!is.finite(m)] <- 0
  } else {
    m <- rep(0, p)
  }

  if (!is.null(S_init)) {
    S <- as.matrix(S_init)
    if (!all(dim(S) == c(p, p))) stop("S_init dimensions do not match ncol(X)")
  } else {
    S <- diag(init_var, p)
  }
  S <- make_pd(S, min_eigen = min_var)
  warm_started <- !is.null(m_init) && !is.null(S_init)

  if (!is.null(tau0_2_init)) {
    tau0_2 <- as.numeric(tau0_2_init)
  } else if (warm_started) {
    agp_var <- pmin(pmax(m^2 + diag(S), min_var), max_var)
    tau0_2 <- min(max(mean(agp_var), min_global_var), max_global_var)
  } else {
    tau0_2 <- 1
  }

  if (!is.null(tau2_init)) {
    tau2 <- as.numeric(tau2_init)
    if (length(tau2) != p) stop("tau2_init length does not match ncol(X)")
  } else if (warm_started) {
    agp_var <- pmin(pmax(m^2 + diag(S), min_var), max_var)
    tau2 <- pmin(pmax(agp_var / sqrt(tau0_2), min_var), max_var)
  } else {
    tau2 <- rep(0.5, p)
  }

  tau0_2 <- min(max(tau0_2, min_global_var), max_global_var)
  tau2 <- pmin(pmax(tau2, min_var), max_var)

  prev_elbo <- gl_blended_igh_elbo(
    X, time, status, m, S, tau2, tau0_2, a, b, a0, b0
  )
  if (!is.finite(prev_elbo)) {
    stop("Initial GL-IGH-VI objective is not finite; check scaling/initialization.")
  }

  elbo_trace <- rep(NA_real_, max_iter)
  direction_norm_trace <- rep(NA_real_, max_iter)
  step_trace <- rep(NA_real_, max_iter)
  n_fallback_used <- 0L
  n_mean_only_rescue_used <- 0L
  n_raw_gradient_rescue_used <- 0L
  update_mode_trace <- rep(NA_character_, max_iter)
  small_change_streak <- 0L
  failed_step_streak <- 0L
  max_failed_step_streak <- 0L
  convergence_status <- "max_iterations"

  for (iter in seq_len(max_iter)) {
    effective_var <- pmax(sqrt(tau0_2) * tau2, min_var)
    gh <- cox_gradient_hessian(X, time, status, m, S)
    gradient <- gh$gradient - m / effective_var

    precision <- gh$hessian +
      diag(1 / effective_var, p)

    S_proposed <- block_solve_precision(
      precision = precision,
      block_size = block_size,
      ridge = ridge,
      min_var = min_var,
      max_var = max_var,
      corr_shrinkage = corr_shrinkage
    )

    # Approximate Newton / natural-gradient direction. S_proposed approximates
    # precision^{-1}, so this preconditions the score by local curvature.
    # Only the global L2 norm is capped; relative coefficient information is
    # preserved, avoiding the plateaus caused by elementwise clipping.
    direction <- as.numeric(S_proposed %*% gradient)
    direction[!is.finite(direction)] <- 0
    dnorm <- sqrt(sum(direction^2))
    if (is.finite(dnorm) && dnorm > max_direction_norm) {
      direction <- direction * (max_direction_norm / dnorm)
      dnorm <- max_direction_norm
    }
    direction_norm_trace[iter] <- dnorm

    alpha <- 1
    accepted <- FALSE
    slack <- 1e-8 * (abs(prev_elbo) + 1)

    for (bt in seq_len(max_backtrack)) {
      step_fraction <- min(1, eta * alpha)
      m_candidate <- m + step_fraction * direction

      # Coherent full-state line search: all candidate quantities approach
      # the current state as step_fraction -> 0.
      S_candidate <- (1 - step_fraction) * S + step_fraction * S_proposed
      S_candidate <- make_pd(S_candidate, min_eigen = min_var)
      second_moment <- m_candidate^2 + diag(S_candidate)

      agp_part <- second_moment
      igh_part <- (
        2 * b + second_moment / pmax(sqrt(tau0_2), sqrt(min_global_var))
      ) / (2 * a + 1)
      tau2_raw <- blend_weight * agp_part + (1 - blend_weight) * igh_part
      tau2_raw <- pmin(pmax(tau2_raw, min_var), max_var)

      tau0_2_raw <- (
        2 * b0 + sum(second_moment / pmax(tau2_raw, min_var))
      ) / (2 * a0 + p)
      tau0_2_raw <- min(max(tau0_2_raw, min_global_var), max_global_var)

      local_damp <- min(1, variance_damping * step_fraction)
      tau2_candidate <- (1 - local_damp) * tau2 + local_damp * tau2_raw
      tau2_candidate <- pmin(pmax(tau2_candidate, min_var), max_var)
      tau0_2_candidate <- (1 - local_damp) * tau0_2 + local_damp * tau0_2_raw
      tau0_2_candidate <- min(max(tau0_2_candidate, min_global_var), max_global_var)

      elbo_candidate <- gl_blended_igh_elbo(
        X, time, status,
        m_candidate, S_candidate,
        tau2_candidate, tau0_2_candidate,
        a, b, a0, b0
      )

      if (is.finite(elbo_candidate) && elbo_candidate >= prev_elbo - slack) {
        m_new <- m_candidate
        S_new <- S_candidate
        tau2_new <- tau2_candidate
        tau0_2_new <- tau0_2_candidate
        elbo <- elbo_candidate
        accepted <- TRUE
        update_mode_trace[iter] <- "full_state"
        step_trace[iter] <- step_fraction
        break
      }
      alpha <- alpha / 2
    }

    # ------------------------------------------------------------------
    # Safeguarded rescue ladder
    #
    # A joint mean/covariance/hypervariance proposal can occasionally fail
    # even when the mean direction itself is an ascent direction.  Rather
    # than declaring a no-op fallback immediately, first try a monotone
    # mean-only line search with S, tau2 and tau0_2 held fixed.  If the
    # block-preconditioned direction is still rejected, try the raw score
    # direction (globally norm-capped, never elementwise clipped).
    #
    # This preserves the GL-IGH objective and avoids the DBCD pathology in
    # which every joint proposal was rejected after an AGP warm start.
    # ------------------------------------------------------------------
    if (!accepted && isTRUE(rescue_mean_only)) {
      rescue_alpha <- 1
      rescue_gain <- min_rescue_gain_rel * (abs(prev_elbo) + 1)

      for (bt in seq_len(rescue_max_backtrack)) {
        rescue_step <- min(1, eta * rescue_alpha)
        m_candidate <- m + rescue_step * direction

        elbo_candidate <- gl_blended_igh_elbo(
          X, time, status,
          m_candidate, S, tau2, tau0_2,
          a, b, a0, b0
        )

        if (is.finite(elbo_candidate) &&
            elbo_candidate >= prev_elbo + rescue_gain) {
          m_new <- m_candidate
          S_new <- S
          tau2_new <- tau2
          tau0_2_new <- tau0_2
          elbo <- elbo_candidate
          accepted <- TRUE
          n_mean_only_rescue_used <- n_mean_only_rescue_used + 1L
          update_mode_trace[iter] <- "mean_only_rescue"
          step_trace[iter] <- rescue_step
          break
        }

        rescue_alpha <- rescue_alpha / 2
      }
    }

    if (!accepted && isTRUE(rescue_raw_gradient)) {
      raw_direction <- as.numeric(gradient)
      raw_direction[!is.finite(raw_direction)] <- 0
      raw_norm <- sqrt(sum(raw_direction^2))

      if (is.finite(raw_norm) && raw_norm > 0) {
        if (raw_norm > max_direction_norm) {
          raw_direction <- raw_direction * (max_direction_norm / raw_norm)
        }

        rescue_alpha <- 1
        rescue_gain <- min_rescue_gain_rel * (abs(prev_elbo) + 1)

        for (bt in seq_len(rescue_max_backtrack)) {
          rescue_step <- min(1, eta * rescue_alpha)
          m_candidate <- m + rescue_step * raw_direction

          elbo_candidate <- gl_blended_igh_elbo(
            X, time, status,
            m_candidate, S, tau2, tau0_2,
            a, b, a0, b0
          )

          if (is.finite(elbo_candidate) &&
              elbo_candidate >= prev_elbo + rescue_gain) {
            m_new <- m_candidate
            S_new <- S
            tau2_new <- tau2
            tau0_2_new <- tau0_2
            elbo <- elbo_candidate
            accepted <- TRUE
            n_raw_gradient_rescue_used <- n_raw_gradient_rescue_used + 1L
            update_mode_trace[iter] <- "raw_gradient_rescue"
            step_trace[iter] <- rescue_step
            break
          }

          rescue_alpha <- rescue_alpha / 2
        }
      }
    }

    if (!accepted) {
      n_fallback_used <- n_fallback_used + 1L
      failed_step_streak <- failed_step_streak + 1L
      max_failed_step_streak <- max(max_failed_step_streak, failed_step_streak)
      m_new <- m
      S_new <- S
      tau2_new <- tau2
      tau0_2_new <- tau0_2
      elbo <- prev_elbo
      update_mode_trace[iter] <- "fallback"
      step_trace[iter] <- 0
    } else {
      failed_step_streak <- 0L
    }

    elbo_trace[iter] <- elbo

    if (iter > 1 && accepted) {
      rel_change <- abs(elbo - prev_elbo) / (abs(prev_elbo) + 1e-8)
      if (is.finite(rel_change) && rel_change < tol) {
        small_change_streak <- small_change_streak + 1L
      } else {
        small_change_streak <- 0L
      }
    } else {
      small_change_streak <- 0L
    }

    m <- m_new
    S <- S_new
    tau2 <- tau2_new
    tau0_2 <- tau0_2_new
    prev_elbo <- elbo

    if (verbose && (iter == 1 || iter %% 25 == 0)) {
      message(
        "GL-IGH-VI iter=", iter,
        ", objective=", round(elbo, 4),
        ", tau0_2=", signif(tau0_2, 4),
        ", ||direction||=", signif(direction_norm_trace[iter], 4),
        ", step=", signif(step_trace[iter], 4)
      )
    }

    if (small_change_streak >= patience) {
      convergence_status <- "converged"
      elbo_trace <- elbo_trace[seq_len(iter)]
      direction_norm_trace <- direction_norm_trace[seq_len(iter)]
      step_trace <- step_trace[seq_len(iter)]
      update_mode_trace <- update_mode_trace[seq_len(iter)]
      break
    }

    if (failed_step_streak >= stall_patience) {
      convergence_status <- "stalled"
      elbo_trace <- elbo_trace[seq_len(iter)]
      direction_norm_trace <- direction_norm_trace[seq_len(iter)]
      step_trace <- step_trace[seq_len(iter)]
      update_mode_trace <- update_mode_trace[seq_len(iter)]
      break
    }
  }

  list(
    method = "Global-Local Gaussian–IGH-VI",
    beta = m,
    covariance = S,
    tau2 = tau2,
    tau0_2 = tau0_2,
    effective_var = sqrt(tau0_2) * tau2,
    a = a,
    b = b,
    a0 = a0,
    b0 = b0,
    blend_weight = blend_weight,
    elbo = elbo_trace[is.finite(elbo_trace)],
    iterations = sum(is.finite(elbo_trace)),
    n_fallback_used = n_fallback_used,
    n_mean_only_rescue_used = n_mean_only_rescue_used,
    n_raw_gradient_rescue_used = n_raw_gradient_rescue_used,
    max_failed_step_streak = max_failed_step_streak,
    convergence_status = convergence_status,
    warm_started = warm_started,
    tol = tol,
    patience = patience,
    stall_patience = stall_patience,
    direction_norm = direction_norm_trace[is.finite(direction_norm_trace)],
    accepted_step = step_trace[is.finite(step_trace)],
    update_mode = update_mode_trace[!is.na(update_mode_trace)],
    optimizer = paste(
      "block-preconditioned variational mean update with full-state backtracking;",
      "safeguarded mean-only and raw-gradient rescue line searches"
    )
  )
}

gl_blended_igh_elbo <- function(X,
                                time,
                                status,
                                m,
                                S,
                                tau2,
                                tau0_2,
                                a,
                                b,
                                a0,
                                b0) {
  
  tau2 <- pmax(tau2, 1e-8)
  
  tau0_2 <- max(tau0_2, 1e-8)
  
  effective_var <- pmax(
    sqrt(tau0_2) * tau2,
    1e-8
  )
  
  # ------------------------------------------------------------
  # Expected Cox partial log-likelihood
  # ------------------------------------------------------------
  
  expected_loglik <- cox_expected_loglik(
    X = X,
    time = time,
    status = status,
    m = m,
    S = S
  )
  
  if (!is.finite(expected_loglik)) {
    return(NA_real_)
  }
  
  # ------------------------------------------------------------
  # Gaussian prior contribution
  # ------------------------------------------------------------
  
  expected_logprior_beta <- -0.5 * sum(
    log(2 * pi * effective_var) +
      (m^2 + diag(S)) / effective_var
  )
  
  # ------------------------------------------------------------
  # Local inverse-gamma prior
  # ------------------------------------------------------------
  
  logprior_tau2 <- sum(
    a * log(b) -
      lgamma(a) -
      (a + 1) * log(tau2) -
      b / tau2
  )
  
  # ------------------------------------------------------------
  # Global inverse-gamma prior
  # ------------------------------------------------------------
  
  logprior_tau0 <- (
    a0 * log(b0) -
      lgamma(a0) -
      (a0 + 1) * log(tau0_2) -
      b0 / tau0_2
  )
  
  # ------------------------------------------------------------
  # Entropy term
  # ------------------------------------------------------------
  
  entropy <- gaussian_entropy(S)
  
  # ------------------------------------------------------------
  # ELBO
  # ------------------------------------------------------------
  
  elbo <- (
    expected_loglik +
      expected_logprior_beta +
      logprior_tau2 +
      logprior_tau0 +
      entropy
  )
  
  if (!is.finite(elbo)) {
    return(NA_real_)
  }
  
  return(elbo)
}

# ============================================================================
# adaptive wrapper for GL-IGH-VI
# ============================================================================

fit_gl_blended_igh_vi_adaptive <- function(
    X,
    time,
    status,
    agp_fit = NULL,
    temper_scales = c(0, 0.25, 0.50, 0.75, 1.00),
    include_cold = TRUE,
    multistart_mode = c("adaptive", "always"),
    fallback_fraction_trigger = 0.10,
    verbose_multistart = TRUE,
    ...) {
  
  multistart_mode <- match.arg(multistart_mode)
  
  X <- as.matrix(X)
  p <- ncol(X)
  
  # ============================================================
  # Helper: safely extract ELBO trace
  # ============================================================
  
  safe_elbo <- function(fit) {
    if (is.null(fit) ||
        is.null(fit$elbo) ||
        length(fit$elbo) == 0) {
      return(numeric(0))
    }
    
    elbo <- as.numeric(fit$elbo)
    elbo[is.finite(elbo)]
  }
  
  # ============================================================
  # Helper: safely extract final ELBO
  # ============================================================
  
  safe_final_elbo <- function(fit) {
    
    elbo <- safe_elbo(fit)
    
    if (length(elbo) == 0) {
      return(NA_real_)
    }
    
    as.numeric(tail(elbo, 1))
  }
  
  # ============================================================
  # Helper: fit one candidate
  # ============================================================
  
  fit_candidate <- function(
    label,
    m_init = NULL,
    S_init = NULL) {
    
    ans <- tryCatch(
      
      fit_gl_blended_igh_vi(
        X = X,
        time = time,
        status = status,
        m_init = m_init,
        S_init = S_init,
        ...
      ),
      
      error = function(e) {
        
        structure(
          list(
            method =
              "Global-Local Gaussian–IGH-VI",
            
            beta =
              rep(NA_real_, p),
            
            covariance =
              diag(NA_real_, p),
            
            elbo =
              numeric(),
            
            iterations =
              0L,
            
            n_fallback_used =
              NA_integer_,
            
            n_mean_only_rescue_used =
              NA_integer_,
            
            n_raw_gradient_rescue_used =
              NA_integer_,
            
            max_failed_step_streak =
              NA_integer_,
            
            convergence_status =
              "error",
            
            error_message =
              conditionMessage(e)
          ),
          
          class = "gl_igh_failed_fit"
        )
      }
    )
    
    ans$start_label <- label
    
    ans
  }
  
  # ============================================================
  # Helper: create diagnostic row for one candidate
  # ============================================================
  
  diagnostic_row <- function(
    fit,
    start_key) {
    
    elbo <- safe_elbo(fit)
    
    n_iter <- length(elbo)
    
    initial_elbo <-
      if (n_iter >= 1) {
        elbo[1]
      } else {
        NA_real_
      }
    
    final_elbo <-
      if (n_iter >= 1) {
        tail(elbo, 1)
      } else {
        NA_real_
      }
    
    elbo_gain <-
      if (n_iter >= 2) {
        final_elbo - initial_elbo
      } else {
        NA_real_
      }
    
    last_rel_change <-
      if (n_iter >= 2) {
        
        abs(
          elbo[n_iter] -
            elbo[n_iter - 1]
        ) /
          (
            abs(elbo[n_iter - 1]) +
              1e-8
          )
        
      } else {
        
        NA_real_
      }
    
    fallbacks <-
      if (!is.null(fit$n_fallback_used)) {
        as.integer(fit$n_fallback_used)
      } else {
        NA_integer_
      }
    
    fallback_fraction <-
      if (is.finite(fallbacks) &&
          n_iter > 0) {
        
        fallbacks / n_iter
        
      } else {
        
        NA_real_
      }
    
    data.frame(
      
      Start_Key =
        start_key,
      
      Start_Label =
        if (!is.null(fit$start_label)) {
          fit$start_label
        } else {
          start_key
        },
      
      Status =
        if (!is.null(fit$convergence_status)) {
          fit$convergence_status
        } else {
          "unknown"
        },
      
      Iterations =
        as.integer(n_iter),
      
      Initial_ELBO =
        as.numeric(initial_elbo),
      
      Final_ELBO =
        as.numeric(final_elbo),
      
      ELBO_Gain =
        as.numeric(elbo_gain),
      
      Last_Rel_Change =
        as.numeric(last_rel_change),
      
      Fallbacks =
        fallbacks,
      
      Fallback_Fraction =
        as.numeric(fallback_fraction),
      
      Mean_Only_Rescues =
        if (!is.null(
          fit$n_mean_only_rescue_used
        )) {
          as.integer(
            fit$n_mean_only_rescue_used
          )
        } else {
          NA_integer_
        },
      
      Raw_Gradient_Rescues =
        if (!is.null(
          fit$n_raw_gradient_rescue_used
        )) {
          as.integer(
            fit$n_raw_gradient_rescue_used
          )
        } else {
          NA_integer_
        },
      
      Max_Failed_Streak =
        if (!is.null(
          fit$max_failed_step_streak
        )) {
          as.integer(
            fit$max_failed_step_streak
          )
        } else {
          NA_integer_
        },
      
      stringsAsFactors = FALSE
    )
  }
  
  # ============================================================
  # 1. Cold start
  # ============================================================
  
  fits <- list()
  
  if (isTRUE(include_cold)) {
    
    fits$cold <-
      fit_candidate(
        label = "cold"
      )
  }
  
  cold_fit <-
    if ("cold" %in% names(fits)) {
      fits$cold
    } else {
      NULL
    }
  
  cold_fallback_fraction <-
    if (!is.null(cold_fit) &&
        !is.null(cold_fit$iterations) &&
        cold_fit$iterations > 0 &&
        !is.null(
          cold_fit$n_fallback_used
        ) &&
        is.finite(
          cold_fit$n_fallback_used
        )) {
      
      cold_fit$n_fallback_used /
        cold_fit$iterations
      
    } else {
      
      NA_real_
    }
  
  cold_is_clean <-
    !is.null(cold_fit) &&
    identical(
      cold_fit$convergence_status,
      "converged"
    ) &&
    (
      is.na(cold_fallback_fraction) ||
        cold_fallback_fraction <=
        fallback_fraction_trigger
    )
  
  # ============================================================
  # 2. Decide whether tempered starts are needed
  # ============================================================
  
  need_tempered <-
    identical(
      multistart_mode,
      "always"
    ) ||
    !cold_is_clean
  
  # ============================================================
  # 3. Tempered AGP starts
  # ============================================================
  
  if (need_tempered) {
    
    if (is.null(agp_fit) ||
        is.null(agp_fit$beta) ||
        is.null(agp_fit$covariance)) {
      
      stop(
        "Adaptive tempered multi-start was triggered, ",
        "but a valid AGP-VI fit with beta and covariance ",
        "was not supplied."
      )
    }
    
    agp_beta <-
      as.numeric(
        agp_fit$beta
      )
    
    agp_cov <-
      as.matrix(
        agp_fit$covariance
      )
    
    if (length(agp_beta) != p ||
        !all(
          dim(agp_cov) == c(p, p)
        )) {
      
      stop(
        "AGP initialization dimensions do not match X."
      )
    }
    
    temper_scales <-
      sort(
        unique(
          as.numeric(
            temper_scales
          )
        )
      )
    
    temper_scales <-
      temper_scales[
        is.finite(temper_scales) &
          temper_scales >= 0 &
          temper_scales <= 1
      ]
    
    if (length(temper_scales) == 0) {
      
      stop(
        "No valid temper_scales in [0, 1]."
      )
    }
    
    for (cc in temper_scales) {
      
      key <-
        paste0(
          "agp_",
          format(
            cc,
            trim = TRUE,
            scientific = FALSE
          )
        )
      
      key <-
        gsub(
          "\\.",
          "_",
          key
        )
      
      label <-
        sprintf(
          "%.2f x AGP mean",
          cc
        )
      
      if (isTRUE(verbose_multistart)) {
        
        cat(
          sprintf(
            "    [GL-IGH multi-start] trying %s\n",
            label
          )
        )
      }
      
      fits[[key]] <-
        fit_candidate(
          label = label,
          m_init = cc * agp_beta,
          S_init = agp_cov
        )
    }
  }
  
  # ============================================================
  # 4. Build diagnostic table
  # ============================================================
  
  diagnostics <-
    do.call(
      rbind,
      lapply(
        names(fits),
        function(nm) {
          
          diagnostic_row(
            fit = fits[[nm]],
            start_key = nm
          )
        }
      )
    )
  
  rownames(diagnostics) <- NULL
  
  # ============================================================
  # Print complete candidate comparison
  # ============================================================
  
  if (isTRUE(verbose_multistart)) {
    
    cat(
      "\n============================================\n"
    )
    
    cat(
      "GL-IGH-VI MULTISTART SUMMARY\n"
    )
    
    cat(
      "============================================\n"
    )
    
    print(
      diagnostics,
      row.names = FALSE,
      digits = 6
    )
    
    cat(
      "============================================\n\n"
    )
  }
  
  # ============================================================
  # 5. Select best genuinely converged fit
  # ============================================================
  
  converged <-
    which(
      diagnostics$Status ==
        "converged" &
        is.finite(
          diagnostics$Final_ELBO
        )
    )
  
  if (length(converged) > 0) {
    
    best_row <-
      converged[
        which.max(
          diagnostics$Final_ELBO[
            converged
          ]
        )
      ]
    
    selection_basis <-
      "highest_final_ELBO_among_converged"
    
  } else {
    
    finite_rows <-
      which(
        is.finite(
          diagnostics$Final_ELBO
        )
      )
    
    if (length(finite_rows) == 0) {
      
      stop(
        "All GL-IGH-VI starts failed ",
        "to produce a finite objective."
      )
    }
    
    best_row <-
      finite_rows[
        which.max(
          diagnostics$Final_ELBO[
            finite_rows
          ]
        )
      ]
    
    selection_basis <-
      "highest_finite_ELBO_no_start_converged"
    
    warning(
      "No GL-IGH-VI initialization converged. ",
      "Returning the highest finite objective fit ",
      "for diagnostics only; do not treat it as converged."
    )
  }
  
  # ============================================================
  # 6. Extract selected fit
  # ============================================================
  
  best_key <-
    diagnostics$Start_Key[
      best_row
    ]
  
  best_fit <-
    fits[[best_key]]
  
  # ============================================================
  # 7. Attach multistart information
  # ============================================================
  
  best_fit$multistart_used <-
    need_tempered
  
  best_fit$multistart_mode <-
    multistart_mode
  
  best_fit$selected_start_key <-
    best_key
  
  best_fit$selected_start_label <-
    diagnostics$Start_Label[
      best_row
    ]
  
  best_fit$multistart_selection_basis <-
    selection_basis
  
  best_fit$multistart_diagnostics <-
    diagnostics
  
  # Also provide intuitive alias
  best_fit$multistart_table <-
    diagnostics
  
  best_fit$adaptive_restart_used <-
    need_tempered
  
  # ============================================================
  # 8. Print selected solution
  # ============================================================
  
  if (isTRUE(verbose_multistart)) {
    
    cat(
      sprintf(
        paste0(
          "    [GL-IGH multi-start] selected %s; ",
          "status=%s; iterations=%d; ",
          "final ELBO=%.6f\n"
        ),
        best_fit$selected_start_label,
        best_fit$convergence_status,
        length(
          safe_elbo(best_fit)
        ),
        safe_final_elbo(
          best_fit
        )
      )
    )
    
    cat(
      sprintf(
        "    [GL-IGH multi-start] selection basis: %s\n",
        selection_basis
      )
    )
  }
  
  best_fit
}