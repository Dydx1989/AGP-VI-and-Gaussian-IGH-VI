load_required_packages <- function() {
  
  # --------------------------------------------------------------------------
  # Standard CRAN packages
  # --------------------------------------------------------------------------
  
  required_cran <- c(
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
  
  missing_cran <- required_cran[
    !sapply(
      required_cran,
      requireNamespace,
      quietly = TRUE
    )
  ]
  
  if (length(missing_cran) > 0) {
    
    message(
      "Installing missing CRAN packages: ",
      paste(missing_cran, collapse = ", ")
    )
    
    install.packages(missing_cran)
  }
  
  
  # --------------------------------------------------------------------------
  # rstanarm
  #
  # Install from Stan R-universe because the Bayesian Lasso survival
  # benchmark requires rstanarm::stan_surv().
  # --------------------------------------------------------------------------
  
  rstanarm_ok <- requireNamespace(
    "rstanarm",
    quietly = TRUE
  )
  
  stan_surv_ok <- FALSE
  
  if (rstanarm_ok) {
    stan_surv_ok <- "stan_surv" %in%
      getNamespaceExports("rstanarm")
  }
  
  if (!rstanarm_ok || !stan_surv_ok) {
    
    message(
      "Installing/updating rstanarm from Stan R-universe ",
      "to obtain stan_surv()..."
    )
    
    install.packages(
      "rstanarm",
      repos = c(
        "https://stan-dev.r-universe.dev",
        getOption("repos")
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Final package checks
  # --------------------------------------------------------------------------
  
  required <- c(
    required_cran,
    "rstanarm"
  )
  
  still_missing <- required[
    !sapply(
      required,
      requireNamespace,
      quietly = TRUE
    )
  ]
  
  if (length(still_missing) > 0) {
    stop(
      "The following packages could not be installed: ",
      paste(still_missing, collapse = ", ")
    )
  }
  
  # Make sure stan_surv() is actually available
  if (!"stan_surv" %in% getNamespaceExports("rstanarm")) {
    stop(
      paste0(
        "rstanarm is installed, but stan_surv() is not available. ",
        "Please install the Stan R-universe version of rstanarm ",
        "and restart R."
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Load packages
  # --------------------------------------------------------------------------
  
  suppressPackageStartupMessages({
    library(survival)
    library(glmnet)
    library(rstanarm)
    library(MASS)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(gridExtra)
    library(Boruta)
    library(randomForest)
  })
  
  
  # --------------------------------------------------------------------------
  # Confirmation
  # --------------------------------------------------------------------------
  
  message(
    "All required packages loaded successfully."
  )
  
  message(
    "rstanarm version: ",
    as.character(packageVersion("rstanarm"))
  )
  
  message(
    "stan_surv available: ",
    "stan_surv" %in% getNamespaceExports("rstanarm")
  )
  
  invisible(TRUE)
}