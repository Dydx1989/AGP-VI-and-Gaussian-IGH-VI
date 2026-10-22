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
                       patience = 5,
                       verbose = FALSE) {
  X <- as.matrix(X)
  p <- ncol(X)
  m <- rep(0, p)
  S <- diag(init_var, p)
  sigma2 <- rep(0.5, p)
  elbo_trace <- rep(NA_real_, max_iter)
  # FIX (see PATCH_NOTES_GL_IGH_VI_oscillation.md): track the last accepted
  # ELBO explicitly and require every accepted step to not decrease it
  # (up to a small numerical slack). The original backtracking loop below
  # only checked is.finite(elbo_candidate), which almost never rejects a
  # step in practice (NaN/Inf is rare) and therefore provided no real
  # protection against overshoot. This restores a genuine sufficient-decrease
  # (Armijo-style) line search: an oscillating or overshooting step is now
  # rejected and the step size is halved (up to 15 times) until an
  # ELBO-improving step is found.
  prev_elbo <- -Inf
  n_fallback_used <- 0
  small_change_streak <- 0
  # FIX, ATTEMPT 3 (2026-08-07; applied here for consistency with the
  # identical fix in fit_gl_blended_igh_vi() -- see
  # PATCH_NOTES_GL_IGH_VI_oscillation.md). This function had the SAME
  # unguarded-fallback bug: when backtracking was exhausted, the code fell
  # back to recomputing S via a different, cruder method and accepted the
  # resulting (m, S, sigma2, elbo) UNCONDITIONALLY, with no check against
  # prev_elbo. AGP-VI's per-coordinate sigma2 = m_j^2 + S_jj feedback makes
  # this fire less often and matter less than it did for the global-variance
  # coupling in the GL-IGH model, but it is not immune, especially at very
  # high p (e.g. VDX, p ~ 1000 after screening) where backtracking exhaustion
  # is routine. An unconditionally-accepted fallback state here can silently
  # corrupt m, potentially explaining anomalously shrunk coefficients (e.g.
  # near-zero everywhere, yielding Mean_Selected = 0 despite non-random
  # discrimination) on the highest-p real datasets. The fix is identical:
  # the fallback is now a true no-op (carry the previous state forward
  # unchanged), and convergence requires `patience` (default 5) consecutive
  # small-relative-change iterations rather than just one.
  for (iter in seq_len(max_iter)) {
    gradient_hessian <- cox_gradient_hessian(X, time, status, m, S)
    gradient <- gradient_hessian$gradient - m / pmax(sigma2, min_var)
    gradient <- pmax(pmin(gradient, grad_clip), -grad_clip)
    precision <- gradient_hessian$hessian + diag(1 / pmax(sigma2, min_var), p) + diag(ridge, p)
    S_new <- block_solve_precision(precision, block_size, ridge, min_var, max_var, corr_shrinkage)
    step <- eta
    accepted <- FALSE
    slack <- 1e-6 * (abs(prev_elbo) + 1)
    if (!is.finite(slack)) slack <- 1e-6
    for (bt in seq_len(15)) {
      m_candidate <- pmax(pmin(m + step * gradient, 10), -10)
      sigma2_candidate <- m_candidate^2 + diag(S_new)
      sigma2_candidate <- pmin(pmax(sigma2_candidate, min_var), max_var)
      elbo_candidate <- agp_elbo(X, time, status, m_candidate, S_new, sigma2_candidate)
      if (is.finite(elbo_candidate) && elbo_candidate >= prev_elbo - slack) {
        m_new <- m_candidate
        sigma2_new <- sigma2_candidate
        elbo <- elbo_candidate
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) {
      # ATTEMPT 3 fix: true no-op, matching fit_gl_blended_igh_vi(). No
      # candidate state is constructed or accepted here -- (m, S, sigma2,
      # elbo) are simply carried forward unchanged from the previous
      # iteration.
      n_fallback_used <- n_fallback_used + 1
      m_new <- m
      S_new <- S
      sigma2_new <- sigma2
      elbo <- prev_elbo
    }
    elbo_trace[iter] <- elbo
    prev_elbo <- elbo
    if (verbose && (iter == 1 || iter %% 25 == 0)) {
      message("Block AGP-VI iter = ", iter, ", ELBO = ", round(elbo, 4))
    }
    if (iter > 2 && is.finite(elbo) && is.finite(elbo_trace[iter - 1])) {
      rel_change <- abs(elbo - elbo_trace[iter - 1]) / (abs(elbo_trace[iter - 1]) + 1e-8)
      # FIX, ATTEMPT 6 (matching fit_gl_blended_igh_vi()): a no-op fallback
      # iteration leaves elbo identical to prev_elbo, so rel_change is
      # trivially 0 and would otherwise be indistinguishable from genuine
      # convergence. A fallback iteration now always resets the streak,
      # regardless of its (trivially small) rel_change.
      if (accepted && is.finite(rel_change) && rel_change < tol) {
        small_change_streak <- small_change_streak + 1
      } else {
        small_change_streak <- 0
      }
      # Patience fix, matching fit_gl_blended_igh_vi(): require several
      # consecutive small-change iterations before declaring convergence.
      if (small_change_streak >= patience) {
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
  list(method = "AGP-VI", beta = m, covariance = S, sigma2 = sigma2, elbo = elbo_trace[is.finite(elbo_trace)], iterations = sum(is.finite(elbo_trace)), n_fallback_used = n_fallback_used)
}
