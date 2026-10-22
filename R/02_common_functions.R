safe_exp <- function(x, cap = 30) {
  exp(pmin(pmax(x, -cap), cap))
}

safe_logsumexp <- function(z) {
  zmax <- max(z, na.rm = TRUE)
  if (!is.finite(zmax)) {
    return(NA_real_)
  }
  zmax + log(sum(exp(z - zmax)) + 1e-12)
}

make_pd <- function(M, min_eigen = 1e-6) {
  M <- 0.5 * (M + t(M))
  eig <- eigen(M, symmetric = TRUE)
  eig$values <- pmax(eig$values, min_eigen)
  out <- eig$vectors %*% diag(eig$values, length(eig$values)) %*% t(eig$vectors)
  0.5 * (out + t(out))
}

block_solve_precision <- function(precision,
                                  block_size = 10,
                                  ridge = 0.1,
                                  min_var = 1e-4,
                                  max_var = 2,
                                  corr_shrinkage = 0.5,
                                  min_eigen = 1e-6) {
  p <- ncol(precision)
  S_new <- matrix(0, p, p)
  blocks <- split(seq_len(p), ceiling(seq_len(p) / block_size))
  for (block in blocks) {
    P_block <- precision[block, block, drop = FALSE]
    P_block <- 0.5 * (P_block + t(P_block))
    P_block <- P_block + diag(ridge, length(block))
    eig <- eigen(P_block, symmetric = TRUE)
    eig$values <- pmax(eig$values, ridge)
    S_block <- eig$vectors %*% diag(1 / eig$values, length(eig$values)) %*% t(eig$vectors)
    diag_values <- diag(S_block)
    S_block <- corr_shrinkage * S_block
    diag(S_block) <- diag_values
    diag(S_block) <- pmin(pmax(diag(S_block), min_var), max_var)
    S_block <- make_pd(S_block, min_eigen = min_eigen)
    diag(S_block) <- pmin(pmax(diag(S_block), min_var), max_var)
    S_new[block, block] <- S_block
  }
  S_new
}

cox_order <- function(X, time, status) {
  order_id <- order(time)
  list(X = X[order_id, , drop = FALSE], time = time[order_id], status = status[order_id])
}

cox_expected_loglik <- function(X, time, status, m, S) {
  ordered <- cox_order(X, time, status)
  X <- ordered$X
  time <- ordered$time
  status <- ordered$status
  eta <- drop(X %*% m)
  quad <- rowSums((X %*% S) * X)
  z <- pmin(pmax(eta + 0.5 * quad, -30), 30)
  event_id <- which(status == 1)
  value <- 0
  for (i in event_id) {
    risk <- which(time >= time[i])
    value <- value + eta[i] - safe_logsumexp(z[risk])
  }
  value
}

cox_partial_loglik <- function(X, time, status, beta) {
  eta <- drop(as.matrix(X) %*% beta)
  eta <- pmin(pmax(eta, -30), 30)
  event_id <- which(status == 1)
  loglik <- 0
  for (i in event_id) {
    risk <- which(time >= time[i])
    loglik <- loglik + eta[i] - safe_logsumexp(eta[risk])
  }
  loglik
}

cox_gradient_hessian <- function(X, time, status, m, S) {
  ordered <- cox_order(X, time, status)
  X <- ordered$X
  time <- ordered$time
  status <- ordered$status
  p <- ncol(X)
  eta <- drop(X %*% m)
  quad <- rowSums((X %*% S) * X)
  z <- pmin(pmax(eta + 0.5 * quad, -30), 30)
  w <- safe_exp(z)
  gradient <- rep(0, p)
  hessian <- matrix(0, p, p)
  event_id <- which(status == 1)
  for (i in event_id) {
    risk <- which(time >= time[i])
    X_risk <- X[risk, , drop = FALSE]
    w_risk <- w[risk]
    denominator <- sum(w_risk) + 1e-12
    x_bar <- drop(crossprod(X_risk, w_risk) / denominator)
    second_moment <- crossprod(X_risk, X_risk * w_risk) / denominator
    covariance_risk <- second_moment - tcrossprod(x_bar)
    gradient <- gradient + X[i, ] - x_bar
    hessian <- hessian + covariance_risk
  }
  gradient[!is.finite(gradient)] <- 0
  hessian[!is.finite(hessian)] <- 0
  list(gradient = gradient, hessian = 0.5 * (hessian + t(hessian)))
}

gaussian_entropy <- function(S) {
  S <- make_pd(S)
  p <- ncol(S)
  chol_S <- chol(S)
  logdet <- 2 * sum(log(diag(chol_S)))
  0.5 * (p * (1 + log(2 * pi)) + logdet)
}

agp_elbo <- function(X, time, status, m, S, sigma2) {
  sigma2 <- pmax(sigma2, 1e-8)
  expected_loglik <- cox_expected_loglik(X, time, status, m, S)
  if (!is.finite(expected_loglik)) {
    return(NA_real_)
  }
  expected_logprior <- -0.5 * sum(log(2 * pi * sigma2) + (m^2 + diag(S)) / sigma2)
  output <- expected_loglik + expected_logprior + gaussian_entropy(S)
  ifelse(is.finite(output), output, NA_real_)
}

igh_elbo <- function(X, time, status, m, S, tau2, a, b) {
  tau2 <- pmax(tau2, 1e-8)
  expected_loglik <- cox_expected_loglik(X, time, status, m, S)
  if (!is.finite(expected_loglik)) {
    return(NA_real_)
  }
  expected_logprior_beta <- -0.5 * sum(log(2 * pi * tau2) + (m^2 + diag(S)) / tau2)
  logprior_tau <- sum(a * log(b) - lgamma(a) - (a + 1) * log(tau2) - b / tau2)
  output <- expected_loglik + expected_logprior_beta + logprior_tau + gaussian_entropy(S)
  ifelse(is.finite(output), output, NA_real_)
}
