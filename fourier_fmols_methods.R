fourier_fmols_manual <- function(y, x, bandwidth, k_candidates = 1:5) {
  y <- as.matrix(y)
  x0 <- as.matrix(x)
  T_total <- nrow(y)
  t_index <- 1:T_total

  # ------------------------------------------------
  # PHASE 1: Find Optimal k using OLS & AIC
  # ------------------------------------------------
  best_aic <- Inf
  best_k <- NA
  best_ols_u <- NULL
  best_X_full <- NULL

  for (k in k_candidates) {
    # Generate Fourier terms
    sin_term <- sin(2 * pi * k * t_index / T_total)
    cos_term <- cos(2 * pi * k * t_index / T_total)

    # X_full order: [1=Intercept, 2=Sin, 3=Cos, 4=X]
    X_full_k <- cbind(1, sin_term, cos_term, x0)

    # OLS Estimation
    beta_ols <- solve(t(X_full_k) %*% X_full_k) %*% t(X_full_k) %*% y
    u_hat <- y - X_full_k %*% beta_ols

    # AIC calculation
    ssr <- sum(u_hat^2)
    p <- ncol(X_full_k)
    current_aic <- T_total * log(ssr / T_total) + 2 * p

    # Update if this k is the best so far
    if (current_aic < best_aic) {
      best_aic <- current_aic
      best_k <- k
      best_ols_u <- u_hat
      best_X_full <- X_full_k
    }
  }

  # ------------------------------------------------
  # PHASE 2: FMOLS Correction with Optimal k
  # ------------------------------------------------
  # Calculate first differences of the stochastic variable
  dx <- diff(x0)

  # Trim data to perfectly align with dx (drop first observation)
  u_hat_trim <- best_ols_u[-1, , drop = FALSE]
  x_trim <- best_X_full[-1, , drop = FALSE]
  y_trim <- y[-1, , drop = FALSE]
  T_trim <- nrow(x_trim)

  # Combine error terms
  w <- cbind(u_hat_trim, dx)

  # Compute Long-Run Covariance Matrices using Bartlett Kernel
  Sigma <- (t(w) %*% w) / T_trim
  Lambda <- matrix(0, nrow = ncol(w), ncol = ncol(w))

  for (j in 1:bandwidth) {
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

  # Construct the Corrections
  y_plus <- y_trim - dx %*% solve(Omega_xx) %*% Omega_0x
  Delta_plus <- Delta_0x - Omega_0x %*% solve(Omega_xx) %*% Delta_xx

  # CRITICAL FIX: The correction is a 4x1 matrix!
  # 0 for Intercept, 0 for Sin, 0 for Cos, and Delta_plus for X
  correction <- matrix(c(0, 0, 0, T_trim * as.numeric(Delta_plus)), nrow = 4, ncol = 1)

  # Final FMOLS calculation
  cross_prod_X <- t(x_trim) %*% x_trim
  cross_prod_XY_plus <- t(x_trim) %*% y_plus

  theta <- solve(cross_prod_X) %*% (cross_prod_XY_plus - correction)
  rownames(theta) <- c("Intercept", "Gamma_1 (Sin)", "Gamma_2 (Cos)", "X_Coef")

  # Return the list with optimal k, AIC, and the full coefficient matrix
  return(list(
    Optimal_k = best_k,
    AIC = best_aic,
    Fourier_FMOLS_Estimates = theta
  ))
}
