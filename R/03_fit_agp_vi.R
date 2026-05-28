fit_agp_vi <- function(X,
                       time,
                       status,
                       max_iter = 300,
                       eta = 0.001,
                       tol = 1e-4,
                       init_var = 0.5,
                       min_var = 1e-4,
                       max_var = 5,
                       ridge = 0.05,
                       grad_clip = 5,
                       block_size = 10,
                       corr_shrinkage = 0.3,
                       verbose = FALSE) {
  X <- as.matrix(X)
  p <- ncol(X)
  m <- rep(0, p)
  S <- diag(init_var, p)
  sigma2 <- rep(0.5, p)
  elbo_trace <- rep(NA_real_, max_iter)
  for (iter in seq_len(max_iter)) {
    gradient_hessian <- cox_gradient_hessian(X, time, status, m, S)
    gradient <- gradient_hessian$gradient - m / pmax(sigma2, min_var)
    gradient <- pmax(pmin(gradient, grad_clip), -grad_clip)
    precision <- gradient_hessian$hessian + diag(1 / pmax(sigma2, min_var), p) + diag(ridge, p)
    S_new <- block_solve_precision(precision, block_size, ridge, min_var, max_var, corr_shrinkage)
    step <- eta
    accepted <- FALSE
    for (bt in seq_len(15)) {
      m_candidate <- pmax(pmin(m + step * gradient, 10), -10)
      sigma2_candidate <- m_candidate^2 + diag(S_new)
      sigma2_candidate <- pmin(pmax(sigma2_candidate, min_var), max_var)
      elbo_candidate <- agp_elbo(X, time, status, m_candidate, S_new, sigma2_candidate)
      if (is.finite(elbo_candidate)) {
        m_new <- m_candidate
        sigma2_new <- sigma2_candidate
        elbo <- elbo_candidate
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) {
      m_new <- m
      S_new <- make_pd(S + diag(ridge, p))
      diag(S_new) <- pmin(pmax(diag(S_new), min_var), max_var)
      sigma2_new <- pmin(pmax(m_new^2 + diag(S_new), min_var), max_var)
      elbo <- agp_elbo(X, time, status, m_new, S_new, sigma2_new)
    }
    elbo_trace[iter] <- elbo
    if (verbose && (iter == 1 || iter %% 25 == 0)) {
      message("Block AGP-VI iter = ", iter, ", ELBO = ", round(elbo, 4))
    }
    if (iter > 2 && is.finite(elbo) && is.finite(elbo_trace[iter - 1])) {
      rel_change <- abs(elbo - elbo_trace[iter - 1]) / (abs(elbo_trace[iter - 1]) + 1e-8)
      if (is.finite(rel_change) && rel_change < tol) {
        m <- m_new
        S <- S_new
        sigma2 <- sigma2_new
        elbo_trace <- elbo_trace[seq_len(iter)]
        break
      }
    }
    m <- m_new
    S <- S_new
    sigma2 <- sigma2_new
  }
  list(method = "AGP-VI", beta = m, covariance = S, sigma2 = sigma2, elbo = elbo_trace[is.finite(elbo_trace)], iterations = sum(is.finite(elbo_trace)))
}
