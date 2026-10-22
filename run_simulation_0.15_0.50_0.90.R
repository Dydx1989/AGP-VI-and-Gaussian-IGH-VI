# ============================================================================
# 09_main_simulation.R
# Main simulation study
# ============================================================================

rm(list = ls())

# ----------------------------------------------------------------------------
# Load packages and source functions
# ----------------------------------------------------------------------------

source("R/00_packages.R")
load_required_packages()

source("R/01_simulate_data.R")
source("R/02_common_functions.R")
source("R/03_fit_agp_vi.R")
source("R/04b_fit_gl_blended_igh_vi.R")
source("R/05_fit_existing_methods.R")
source("R/06_evaluation.R")
source("R/07_summary_latex.R")
source("R/08_plots.R")


Sys.setenv(SIM_RHO = "0.15")
source("R/09_main_simulation.R")



Sys.setenv(SIM_RHO = "0.50")
source("R/09_main_simulation.R")



Sys.setenv(SIM_RHO = "0.90")
source("R/09_main_simulation.R")






