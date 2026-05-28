source("R/00_packages.R")
load_required_packages()

source("R/01_simulate_data.R")
source("R/02_common_functions.R")
source("R/03_fit_agp_vi.R")
source("R/04b_fit_gl_blended_igh_vi.R")
source("R/05_fit_existing_methods.R")
source("R/06_evaluation.R")

set.seed(1)

dat <- simulate_cox_data(
  n = 100,
  p = 200,
  n_signal = 10,
  rho = 0.30,
  beta_value = 0.7,
  censoring_rate = 0.50,
  seed = 1
)
split <- train_test_split(
  dat$X,
  dat$time,
  dat$status,
  train_fraction = 0.70,
  seed = 1
)

eval_time <- as.numeric(
  quantile(
    split$time_train[split$status_train == 1],
    0.50
  )
)

fit_cox_lasso <- fit_cox_lasso(
  split$X_train,
  split$time_train,
  split$status_train
)

fit_elasticnet <- fit_cox_elasticnet(
  split$X_train,
  split$time_train,
  split$status_train,
  alpha = 0.5
)

fit_bayes_lasso <- fit_bayesian_lasso_cox(
  split$X_train,
  split$time_train,
  split$status_train,
  chains = 2,
  iter = 500,
  warmup = 250,
  seed = 1
)

fit_agp <- fit_agp_vi(
  split$X_train,
  split$time_train,
  split$status_train,
  max_iter = 50,
  block_size = 10
)

fit_gl_igh <- fit_gl_blended_igh_vi(
  split$X_train,
  split$time_train,
  split$status_train,
  max_iter = 50,
  block_size = 10
)

res_cox_lasso <- evaluate_model(
  fit_cox_lasso,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_elasticnet <- evaluate_model(
  fit_elasticnet,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_bayes_lasso <- evaluate_model(
  fit_bayes_lasso,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_agp <- evaluate_model(
  fit_agp,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_gl_igh <- evaluate_model(
  fit_gl_igh,
  split$X_train,
  split$time_train,
  split$status_train,
  split$X_test,
  split$time_test,
  split$status_test,
  dat$beta_true,
  eval_time = eval_time
)

res_cox_lasso$Method <- "Cox-Lasso"
res_elasticnet$Method <- "Cox-ElasticNet"
res_bayes_lasso$Method <- "Bayesian-Lasso-Cox"
res_agp$Method <- "AGP-VI"
res_gl_igh$Method <- "Global-Local Gaussian–IGH-VI"

quick_results <- rbind(
  res_cox_lasso,
  res_elasticnet,
  res_bayes_lasso,
  res_agp,
  res_gl_igh
)

print(quick_results)
