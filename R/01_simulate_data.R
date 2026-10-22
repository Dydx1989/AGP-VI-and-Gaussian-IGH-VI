simulate_cox_data <- function(n = 100,
                              p = 100,
                              n_signal = 10,
                              rho = 0.15,
                              beta_value = 1.2,
                              censoring_rate = 0.30,
                              seed = 2026) {
  set.seed(seed)
  Sigma <- rho ^ abs(outer(seq_len(p), seq_len(p), "-"))
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  X <- scale(X)
  beta_true <- rep(0, p)
  signal_index <- seq_len(n_signal)
  signal_sign <- rep(c(1, -1), length.out = n_signal)
  beta_true[signal_index] <- signal_sign * beta_value
  eta <- drop(X %*% beta_true)
  eta <- pmin(pmax(eta, -10), 10)
  baseline_rate <- 0.05
  event_time <- rexp(n, rate = baseline_rate * exp(eta))
  censor_rate <- tryCatch(
    uniroot(
      function(cr) {
        mean(rexp(3000, rate = cr) < sample(event_time, 3000, replace = TRUE)) - censoring_rate
      },
      lower = 1e-5,
      upper = 10
    )$root,
    error = function(e) 0.05
  )
  censor_time <- rexp(n, rate = censor_rate)
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)
  list(
    X = as.matrix(X),
    time = time,
    status = status,
    beta_true = beta_true,
    signal_index = signal_index,
    observed_censoring = mean(status == 0)
  )
}

train_test_split <- function(X,
                             time,
                             status,
                             train_fraction = 0.70,
                             seed = 2026) {
  set.seed(seed)
  n <- nrow(X)
  train_id <- sample(seq_len(n), size = floor(train_fraction * n))
  test_id <- setdiff(seq_len(n), train_id)
  X_train <- X[train_id, , drop = FALSE]
  X_test <- X[test_id, , drop = FALSE]
  center <- colMeans(X_train)
  scale_value <- apply(X_train, 2, sd)
  scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
  X_train <- sweep(sweep(X_train, 2, center, "-"), 2, scale_value, "/")
  X_test <- sweep(sweep(X_test, 2, center, "-"), 2, scale_value, "/")
  list(
    X_train = as.matrix(X_train),
    X_test = as.matrix(X_test),
    time_train = time[train_id],
    time_test = time[test_id],
    status_train = status[train_id],
    status_test = status[test_id]
  )
}
