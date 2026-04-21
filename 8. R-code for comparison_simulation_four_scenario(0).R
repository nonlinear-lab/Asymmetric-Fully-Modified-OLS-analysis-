# Setup
setwd("C:/Users/User/Desktop/Fumitaka/R")

# Load all method files
source("dols_methods.R")
source("fourier_dols_methods.R")
source("fmols_methods.R")
source("fourier_fmols_methods.R")
source("asymmetric_dols_methods.R")
source("asymmetric_fmols_methods.R")
source("simulation_methods.R")

# Simulation Parameters
T_val     <- 200
reps      <- 1000
true_beta <- 2.5
bandwidth <- 8
n_leads   <- 2
n_lags    <- 2

scenario_names <- c("Strong (iid)",
                    "Weak (AR=0.8)",
                    "Very Weak (AR=0.95)",
                    "No Cointegration (RW)")

set.seed(42)
for (j in 1:4) {
  cat("\n===", scenario_names[j], "===\n")

  results <- simulation_analysis(
    scenario  = as.character(j),
    n_leads   = n_leads,
    n_lags    = n_lags,
    bandwidth = bandwidth,
    T_val     = T_val,
    reps      = reps,
    true_beta = true_beta
  )

  print(round(results, 4))              # ✅ print correct object
}
