fmols_manual <- function(y, x, bandwidth) {
  y <- as.matrix(y)
  x0 <- as.matrix(x)

  # 1. Full X matrix with Intercept (D matrix in cointReg)
  x <- cbind(1, x0)

  # 2. Initial OLS on the FULL data
  beta_ols <- solve(t(x) %*% x) %*% t(x) %*% y
  u_hat <- y - x %*% beta_ols

  # 3. Calculate first differences of the stochastic variable
  dx <- diff(x0)

  # 4. Trim data to perfectly align with dx
  u_hat_trim <- u_hat[-1, , drop = FALSE]
  x_trim <- x[-1, , drop = FALSE]
  y_trim <- y[-1, , drop = FALSE]
  T_trim <- nrow(x_trim)

  # 5. Combine error terms (CRITICAL FIX: Do NOT demean this!)
  w <- cbind(u_hat_trim, dx)

  # 6. Compute Long-Run Covariance Matrices using Bartlett Kernel
  Sigma <- (t(w) %*% w) / T_trim
  Lambda <- matrix(0, nrow = ncol(w), ncol = ncol(w))

  for (j in 1:bandwidth) {
    # Standard Bartlett weight without the +1 denominator
    weight <- 1 - (j / bandwidth)
    Gamma_j <- (t(w[(j+1):T_trim, ]) %*% w[1:(T_trim-j), ]) / T_trim
    Lambda <- Lambda + weight * Gamma_j
  }

  Omega <- Sigma + Lambda + t(Lambda)
  Delta <- Sigma + Lambda

  Omega_00 <- Omega[1, 1]
  Omega_0x <- Omega[1, 2]
  Omega_xx <- Omega[2, 2]
  Delta_0x <- Delta[1, 2]
  Delta_xx <- Delta[2, 2]

  # 7. Construct the Corrections
  y_plus <- y_trim - dx %*% solve(Omega_xx) %*% Omega_0x
  Delta_plus <- Delta_0x - Omega_0x %*% solve(Omega_xx) %*% Delta_xx

  # The correction is ONLY applied to the stochastic beta, not delta
  correction <- matrix(c(0, T_trim * as.numeric(Delta_plus)), nrow = 2, ncol = 1)

  # 8. Final FMOLS calculation
  XX <- t(x_trim) %*% x_trim
  Yplus_X <- t(x_trim) %*% y_plus

  theta <- solve(XX) %*% (Yplus_X - correction)

  # Return the stochastic slope (beta)
  return(theta[2, 1])
}