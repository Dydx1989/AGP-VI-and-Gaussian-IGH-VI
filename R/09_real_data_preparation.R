prepare_model_matrix <- function(data, time_col, status_col, event_code = 1, remove_cols = NULL) {
  data <- as.data.frame(data)
  if (!is.null(remove_cols)) data <- data[, !(names(data) %in% remove_cols), drop = FALSE]
  data <- data[complete.cases(data), , drop = FALSE]
  time <- as.numeric(data[[time_col]])
  status <- as.integer(data[[status_col]] == event_code)
  predictors <- data[, !(names(data) %in% c(time_col, status_col)), drop = FALSE]
  predictors <- predictors[, sapply(predictors, function(x) length(unique(x)) > 1), drop = FALSE]
  X <- model.matrix(~ ., data = predictors)[, -1, drop = FALSE]
  keep <- apply(X, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  list(X = as.matrix(X), time = time, status = status)
}

split_real_data <- function(X, time, status, train_fraction = 0.70, seed = 2026) {
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
  list(X_train = as.matrix(X_train), X_test = as.matrix(X_test), time_train = time[train_id], time_test = time[test_id], status_train = status[train_id], status_test = status[test_id])
}

load_real_datasets <- function() {
  data("lung", package = "survival")
  data("pbc", package = "survival")
  data("colon", package = "survival")
  lung_data <- survival::lung
  lung_data$status_event <- ifelse(lung_data$status == 2, 1, 0)
  lung_prepared <- prepare_model_matrix(lung_data, time_col = "time", status_col = "status_event", event_code = 1, remove_cols = c("inst", "status"))
  pbc_data <- survival::pbc
  pbc_data$status_event <- ifelse(pbc_data$status == 2, 1, 0)
  pbc_prepared <- prepare_model_matrix(pbc_data, time_col = "time", status_col = "status_event", event_code = 1, remove_cols = c("id", "status"))
  colon_data <- survival::colon
  colon_data <- colon_data[colon_data$etype == 2, , drop = FALSE]
  colon_prepared <- prepare_model_matrix(colon_data, time_col = "time", status_col = "status", event_code = 1, remove_cols = c("id", "study", "etype"))
  list(Lung = lung_prepared, PBC = pbc_prepared, Colon = colon_prepared)
}
