#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: Asymmetric-FMOLS Monte Carlo (1000 Reps)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/R/fmols")
library(cointReg)

# --- Simulation Parameters ---
T_val <- 200     # Sample size
reps <- 1000     # Number of repetitions
true_beta <- 2.5
bandwidth <- 2   # Fixed bandwidth for comparison

# Matrices to store t-statistics and slopes for 4 scenarios
sim_plus          <- matrix(NA, nrow = reps, ncol = 4)
sim_minus         <- matrix(NA, nrow = reps, ncol = 4)
sim_tstats_plus   <- matrix(NA, nrow = reps, ncol = 4)
sim_tstats_minus  <- matrix(NA, nrow = reps, ncol = 4)
colnames(sim_tstats_plus)  <- c("Strong", "Weak", "Very_Weak", "None")
colnames(sim_tstats_minus) <- c("Strong", "Weak", "Very_Weak", "None")

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: FMLOS function
#++++++++++++++++++++++++++++++++++++++++++++++++
manual_fmols <- function(y, x, bandwidth) {
  y <- as.matrix(y)
  x0 <- as.matrix(x)

  # 1. Full X matrix with Intercept (D matrix in cointReg)
  x <- cbind(1, x0)

  # 2. Calculate first differences of the stochastic variable
  dx <- diff(x0)
  dx_pos <- ifelse(dx > 0, dx, 0)
  dx_neg <- ifelse(dx < 0, dx, 0)
  x_pos <- cumsum(dx_pos)
  x_neg <- cumsum(dx_neg)
  x_asym <- cbind(1, x_pos, x_neg)

  # 3. Trim data to perfectly align with dx
  # u_hat_trim     <- u_hat[-1, , drop = FALSE]
  x_trim         <- cbind(1, x_pos, x_neg)
  y_trim         <- y[-1, , drop = FALSE]
  T_trim         <- nrow(x_trim)

  # 4. Initial OLS on the FULL data
  #beta_ols   <- solve(t(x) %*% x) %*% t(x) %*% y
  #u_hat      <- y - x %*% beta_ols
  beta_ols   <- solve(t(x_trim) %*% x_trim) %*% t(x_trim) %*% y_trim
  u_hat_trim <- y_trim - x_trim %*% beta_ols

  # 5. Combine error terms (CRITICAL FIX: Do NOT demean this!)
  T1  <- nrow(y_trim)
  w <- cbind(u_hat_trim, dx_pos, dx_neg)

  # 6. Compute Long-Run Covariance Matrices using Bartlett Kernel
  # Sigma  <- (t(w) %*% w) / T1
  Sigma  <- (t(w) %*% w) / (T1 - ncol(x_trim))   # T1 - 3
  # Lambda <- matrix(0, nrow = ncol(w), ncol = ncol(w))
  Lambda <- matrix(0, 3, 3)

  for (j in 1:bandwidth) {
    # Standard Bartlett weight without the +1 denominator
    # weight  <- 1 - (j / (bandwidth + 1))
    weight  <- 1 - (j / bandwidth)
    Gamma_j <- (t(w[(j+1):T1, ]) %*% w[1:(T1-j), ]) / T1
    Lambda  <- Lambda + weight * Gamma_j
  }

  Omega <- Sigma + Lambda + t(Lambda)
  Delta <- Sigma + Lambda

  Omega_0x <- Omega[1, 2:3, drop = FALSE]   # 1×2
  Omega_xx <- Omega[2:3, 2:3]               # 2×2
  Delta_0x <- Delta[1, 2:3, drop = FALSE]   # 1×2
  Delta_xx <- Delta[2:3, 2:3]               # 2×2

  # 7. Construct the Corrections
  dx_mat  <- cbind(dx_pos, dx_neg)          # T1 × 2
  # y_plus <- y_trim - dx_mat %*% solve(Omega_xx) %*% Omega_0x
  y_plus <- y_trim - dx_mat %*% solve(Omega_xx) %*% t(Omega_0x)
  Delta_plus <- Delta_0x - Omega_0x %*% solve(Omega_xx) %*% Delta_xx

  # The correction is ONLY applied to the stochastic beta, not delta
  # correction <- matrix(c(0, T_trim * as.numeric(Delta_plus)), nrow = 2, ncol = 1)
  correction <- matrix(c(0, T1*Delta_plus[1,1], T1*Delta_plus[1,2]), nrow = 3, ncol = 1) # 3×1

  # 8. Final FMOLS calculation
  XX <- t(x_trim) %*% x_trim
  Yplus_X <- t(x_trim) %*% y_plus
  XtX_inv <- tryCatch(solve(t(x_trim) %*% x_trim), error = function(e) NULL)

  theta <- solve(XX) %*% (Yplus_X - correction)

  # Calculate Residuals and Conventional Variance
  u_hat <- y_trim - x_trim %*% theta
  df <- nrow(x_trim) - ncol(x_trim)
  sigma_sq <- sum(u_hat^2) / df

  cov_matrix <- sigma_sq * XtX_inv
  se_hat <- sqrt(diag(cov_matrix))

  # Return the stochastic slope (beta)
  results    <- list(
    beta_pos = theta[2, 1],
    beta_neg = theta[3, 1],
    se_plus    = se_hat[2],
    se_minus   = se_hat[3]
  )
  return(results)
}
#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: The Simulation Loop
#++++++++++++++++++++++++++++++++++++++++++++++++
set.seed(42)
cat("Starting 1000 Repetitions...\n")

for (i in 1:reps) {
  # Generate X (Random Walk)
  x <- cumsum(rnorm(T_val))

  # Generate 4 Scenarios
  y1 <- 10 + true_beta * x + rnorm(T_val)                            # Strong
  y2 <- 10 + true_beta * x + arima.sim(list(ar=0.8), T_val)          # Weak
  y3 <- 10 + true_beta * x + arima.sim(list(ar=0.95), T_val)         # Very Weak
  y4 <- 10 + true_beta * x + cumsum(rnorm(T_val))                    # None

  for (j in 1:4) {
    y_current <- get(paste0("y", j))
    res <- manual_fmols(y_current, x, bandwidth)
    sim_plus[i, j]  <- res$beta_pos
    sim_minus[i, j] <- res$beta_neg
    # Calculate and Store T-statistics: t = (beta_hat - true_beta) / SE
    sim_tstats_plus[i, j]  <- (res$beta_pos - true_beta) / res$se_plus
    sim_tstats_minus[i, j] <- (res$beta_neg - true_beta) / res$se_minus
  }
  if (i %% 100 == 0) cat("Repetition:", i, "\n")
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Analysis of Results
#++++++++++++++++++++++++++++++++++++++++++++++++
# Rejection Rate at 5% (Testing H0: Beta = 2.5)
# If H0 is TRUE, this should be 0.05.
power_size1 <- apply(sim_tstats_plus, 2, function(x) mean(abs(x) > 1.96))
power_size2 <- apply(sim_tstats_minus, 2, function(x) mean(abs(x) > 1.96))

cat("\nEmpirical Rejection Rates (H0: Beta = 2.5):\n")
print(power_size1)
print(power_size2)

# Summary of Slopes
cat("\nSummary of Slope Estimates:\n")
print(round(apply(sim_plus,  2, summary), 4))
print(round(apply(sim_minus, 2, summary), 4))

