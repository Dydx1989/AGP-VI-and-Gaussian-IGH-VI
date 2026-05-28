# ============================================================================
# Improved Tempered GL-Blended Gaussian--IGH-VI
# Better prediction + controlled false positives
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
                                  verbose = FALSE) {
  
  X <- as.matrix(X)
  p <- ncol(X)
  
  blend_weight <- min(max(blend_weight, 0), 1)
  
  m <- rep(0, p)
  S <- diag(init_var, p)
  
  tau2 <- rep(0.5, p)
  tau0_2 <- 1
  
  elbo_trace <- rep(NA_real_, max_iter)
  
  for (iter in seq_len(max_iter)) {
    
    effective_var <- pmax(
      sqrt(tau0_2) * tau2,
      min_var
    )
    
    gradient_hessian <- cox_gradient_hessian(
      X,
      time,
      status,
      m,
      S
    )
    
    gradient <- gradient_hessian$gradient -
      m / effective_var
    
    gradient <- pmax(
      pmin(gradient, grad_clip),
      -grad_clip
    )
    
    precision <- gradient_hessian$hessian +
      diag(1 / effective_var, p) +
      diag(ridge, p)
    
    S_new <- block_solve_precision(
      precision = precision,
      block_size = block_size,
      ridge = ridge,
      min_var = min_var,
      max_var = max_var,
      corr_shrinkage = corr_shrinkage
    )
    
    step <- eta
    accepted <- FALSE
    
    for (bt in seq_len(15)) {
      
      m_candidate <- pmax(
        pmin(m + step * gradient, 10),
        -10
      )
      
      second_moment <- m_candidate^2 + diag(S_new)
      
      agp_part <- second_moment
      
      igh_part <- (
        2 * b +
          second_moment / pmax(sqrt(tau0_2), sqrt(min_global_var))
      ) / (2 * a + 1)
      
      tau2_candidate <- blend_weight * agp_part +
        (1 - blend_weight) * igh_part
      
      tau2_candidate <- pmin(
        pmax(tau2_candidate, min_var),
        max_var
      )
      
      tau0_2_candidate <- (
        2 * b0 +
          sum(second_moment / pmax(tau2_candidate, min_var))
      ) / (2 * a0 + p)
      
      tau0_2_candidate <- min(
        max(tau0_2_candidate, min_global_var),
        max_global_var
      )
      
      elbo_candidate <- gl_blended_igh_elbo(
        X = X,
        time = time,
        status = status,
        m = m_candidate,
        S = S_new,
        tau2 = tau2_candidate,
        tau0_2 = tau0_2_candidate,
        a = a,
        b = b,
        a0 = a0,
        b0 = b0
      )
      
      if (is.finite(elbo_candidate)) {
        m_new <- m_candidate
        tau2_new <- tau2_candidate
        tau0_2_new <- tau0_2_candidate
        elbo <- elbo_candidate
        accepted <- TRUE
        break
      }
      
      step <- step / 2
    }
    
    if (!accepted) {
      
      m_new <- m
      
      S_new <- make_pd(
        S + diag(ridge, p)
      )
      
      diag(S_new) <- pmin(
        pmax(diag(S_new), min_var),
        max_var
      )
      
      second_moment <- m_new^2 + diag(S_new)
      
      agp_part <- second_moment
      
      igh_part <- (
        2 * b +
          second_moment / pmax(sqrt(tau0_2), sqrt(min_global_var))
      ) / (2 * a + 1)
      
      tau2_new <- blend_weight * agp_part +
        (1 - blend_weight) * igh_part
      
      tau2_new <- pmin(
        pmax(tau2_new, min_var),
        max_var
      )
      
      tau0_2_new <- (
        2 * b0 +
          sum(second_moment / pmax(tau2_new, min_var))
      ) / (2 * a0 + p)
      
      tau0_2_new <- min(
        max(tau0_2_new, min_global_var),
        max_global_var
      )
      
      elbo <- gl_blended_igh_elbo(
        X = X,
        time = time,
        status = status,
        m = m_new,
        S = S_new,
        tau2 = tau2_new,
        tau0_2 = tau0_2_new,
        a = a,
        b = b,
        a0 = a0,
        b0 = b0
      )
    }
    
    elbo_trace[iter] <- elbo
    
    if (verbose && (iter == 1 || iter %% 25 == 0)) {
      message(
        "Global-Local Gaussian–IGH-VI iter = ",
        iter,
        ", ELBO = ",
        round(elbo, 4),
        ", tau0_2 = ",
        round(tau0_2_new, 4),
        ", mean_tau2 = ",
        round(mean(tau2_new), 4),
        ", selected-like = ",
        sum(abs(m_new) > 0.15)
      )
    }
    
    if (
      iter > 2 &&
      is.finite(elbo) &&
      is.finite(elbo_trace[iter - 1])
    ) {
      rel_change <- abs(elbo - elbo_trace[iter - 1]) /
        (abs(elbo_trace[iter - 1]) + 1e-8)
      
      if (is.finite(rel_change) && rel_change < tol) {
        m <- m_new
        S <- S_new
        tau2 <- tau2_new
        tau0_2 <- tau0_2_new
        elbo_trace <- elbo_trace[seq_len(iter)]
        break
      }
    }
    
    m <- m_new
    S <- S_new
    tau2 <- tau2_new
    tau0_2 <- tau0_2_new
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
    iterations = sum(is.finite(elbo_trace))
  )
}

# ============================================================================
# ELBO for Improved Tempered GL-Blended Gaussian--IGH-VI
# ============================================================================

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