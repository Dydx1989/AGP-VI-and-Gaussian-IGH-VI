load_required_packages <- function() {
  
  required <- c(
    "survival",
    "glmnet",
    "MASS",
    "ggplot2",
    "dplyr",
    "tidyr",
    "gridExtra",
    "Boruta",
    "randomForest"
  )
  
  missing <- required[
    !sapply(
      required,
      requireNamespace,
      quietly = TRUE
    )
  ]
  
  if (length(missing) > 0) {
    install.packages(missing)
  }
  
  library(survival)
  library(glmnet)
  library(MASS)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(Boruta)
  library(randomForest)
  
  invisible(TRUE)
}