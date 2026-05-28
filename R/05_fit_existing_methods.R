# ============================================================================
# Existing benchmark Cox methods
# Cox-Lasso, Cox-ElasticNet, Bayesian Lasso Cox
# ============================================================================

# ----------------------------------------------------------------------------
# Cox-Lasso
# ----------------------------------------------------------------------------

fit_cox_lasso <- function(X, time, status) {
  
  y <- survival::Surv(time, status)
  
  fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = y,
    family = "cox",
    alpha = 1,
    nfolds = 5,
    maxit = 1e6
  )
  
  beta_hat <- as.numeric(
    coef(fit, s = "lambda.min")
  )
  
  list(
    method = "Cox-Lasso",
    beta = beta_hat,
    fit = fit
  )
}


# ----------------------------------------------------------------------------
# Cox-ElasticNet
# ----------------------------------------------------------------------------

fit_cox_elasticnet <- function(X,
                               time,
                               status,
                               alpha = 0.5) {
  
  y <- survival::Surv(time, status)
  
  fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = y,
    family = "cox",
    alpha = alpha,
    nfolds = 5,
    maxit = 1e6
  )
  
  beta_hat <- as.numeric(
    coef(fit, s = "lambda.min")
  )
  
  list(
    method = "Cox-ElasticNet",
    beta = beta_hat,
    fit = fit,
    alpha = alpha
  )
}


# ----------------------------------------------------------------------------
# Bayesian Lasso Cox
# ----------------------------------------------------------------------------

fit_bayesian_lasso_cox <- function(X,
                                   time,
                                   status,
                                   chains = 2,
                                   iter = 500,
                                   warmup = 250,
                                   seed = 2026) {
  
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    
    warning(
      "rstanarm is not installed. Using Cox-Lasso fallback as Bayesian Lasso Cox approximation."
    )
    
    fallback_fit <- fit_cox_lasso(
      X = X,
      time = time,
      status = status
    )
    
    fallback_fit$method <- "Bayesian-Lasso-Cox"
    fallback_fit$fallback <- TRUE
    
    return(fallback_fit)
  }
  
  data_stan <- as.data.frame(X)
  colnames(data_stan) <- paste0("X", seq_len(ncol(X)))
  data_stan$time <- time
  data_stan$status <- status
  
  formula_stan <- as.formula(
    paste0(
      "survival::Surv(time, status) ~ ",
      paste(
        colnames(data_stan)[seq_len(ncol(X))],
        collapse = " + "
      )
    )
  )
  
  fit <- tryCatch(
    suppressWarnings(
      rstanarm::stan_surv(
        formula = formula_stan,
        data = data_stan,
        basehaz = "ms",
        prior = rstanarm::laplace(
          location = 0,
          scale = 1,
          autoscale = TRUE
        ),
        chains = chains,
        iter = iter,
        warmup = warmup,
        seed = seed,
        refresh = 0
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    
    warning(
      "Bayesian Lasso Cox failed. Using Cox-Lasso fallback."
    )
    
    fallback_fit <- fit_cox_lasso(
      X = X,
      time = time,
      status = status
    )
    
    fallback_fit$method <- "Bayesian-Lasso-Cox"
    fallback_fit$fallback <- TRUE
    
    return(fallback_fit)
  }
  
  posterior_draws <- as.matrix(fit)
  
  beta_cols <- grep(
    "^X[0-9]+$",
    colnames(posterior_draws),
    value = TRUE
  )
  
  beta_hat <- rep(0, ncol(X))
  
  if (length(beta_cols) >= ncol(X)) {
    beta_hat <- colMeans(
      posterior_draws[
        ,
        beta_cols[seq_len(ncol(X))],
        drop = FALSE
      ]
    )
  }
  
  beta_hat[!is.finite(beta_hat)] <- 0
  
  list(
    method = "Bayesian-Lasso-Cox",
    beta = beta_hat,
    fit = fit,
    fallback = FALSE,
    prior = "Laplace"
  )
}