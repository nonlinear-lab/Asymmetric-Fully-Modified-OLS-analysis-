# Setup
setwd("C:/Users/User/Desktop/Fumitaka/R")
library(cointReg)

# Load both engines
source("dols_methods.R")
source("fourier_dols_methods.R")
source("fmols_methods.R")
source("fourier_fmols_methods.R")
source("asymmetric_dols_methods.R")
source("asymmetric_fmols_methods.R")

# Simulation Parameters
T_val     <- 200
reps      <- 1000
true_beta <- 2.5
bandwidth <- 8
n_leads   <- 2
n_lags    <- 2

# Storage for results (Slopes for Scenario 3: AR=0.95)
results_slope <- matrix(NA, nrow = reps, ncol = 9)
colnames(results_slope) <- c("OLS", "DOLS", "Fourier DOLS","FMOLS", "Fourier_FMOLS", "Asymmetricr_DOLS_plus", "Asymmetricr_DOLS_minus", "Asymmetricr_FMOLS_plus", "Asymmetricr_FMOLS_minus")

set.seed(42)
for (i in 1:reps) {
  # Data Generation (Very Weak Cointegration / High Persistence)
  x <- cumsum(rnorm(T_val))
  u <- arima.sim(list(ar=0.95), T_val)
  y <- 10 + true_beta * x + u

  # 1. OLS
  results_slope[i, 1] <- coef(lm(y ~ x))[2]

  # 2. DOLS (from dols_method.R)
  results_slope[i, 2] <- dols_manual(y, x, n_leads, n_lags)

  # 3. DOLS (from dols_method.R)
  results_slope[i, 3] <- fourier_dols_manual(y, x, n_leads, n_lags)$Fourier_DOLS_Estimates["X_Coef", 1]

  # 4. FMOLS (from fmols_methods.R)
  results_slope[i, 4] <- fmols_manual(y, x, bandwidth)

  # 5. Fourier-FMOLS (from fourier_fmols_methods.R)
  results_slope[i, 5] <- fourier_fmols_manual(y, x, bandwidth)$Fourier_FMOLS_Estimates["X_Coef", 1]

  # 6. Asymmetric DOLS (from asymmetric_dols_method.R)
  asym_res1 <- asymmetric_dols_manual(y, x, n_leads, n_lags)
  results_slope[i, 6] <- asym_res1$beta_plus
  results_slope[i, 7] <- asym_res1$beta_minus

  # 6. Asymmetric FMOLS (from asymmetric_fmols_method.R)
  asym_res2 <- asymmetric_fmols_manual(y, x, bandwidth)
  results_slope[i, 8] <- asym_res2$beta_pos
  results_slope[i, 9] <- asym_res2$beta_neg


  if (i %% 100 == 0) cat("Repetition", i, "Complete\n")
}

# Final Comparison Table
summary_stats <- t(apply(results_slope, 2, function(x) {
  c(Mean = mean(x), Bias = mean(x) - true_beta, SD = sd(x), RMSE = sqrt(mean((x - true_beta)^2)))
}))

print(round(summary_stats, 4))
