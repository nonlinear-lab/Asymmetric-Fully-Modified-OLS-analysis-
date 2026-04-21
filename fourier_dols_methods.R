fourier_dols_manual <- function(y, x, n_leads, n_lags, k_candidates = 1:5) { 
  y <- as.matrix(y)
  x <- as.matrix(x)
  N_total <- nrow(x)
  dx <- diff(x)

# No.2 Trim dataset to account for leads and lags
  valid_idx <- (n_lags + 2):(N_total - n_leads)
  y_trim <- y[valid_idx, , drop = FALSE]
  x_trim <- x[valid_idx, , drop = FALSE]
  t_index <- valid_idx
  N_trim <- nrow(y_trim) # Actual number of observations used in regression

# No.3. lead and lag of first difference D(L)delta X
  dx_matrix <- NULL
  for (j in (-n_leads):n_lags) {
    dx_shifted <- dx[(valid_idx - 1 - j), , drop = FALSE]
    dx_matrix <- cbind(dx_matrix, dx_shifted)
  }

# No4. Using AIC to choose optimal frequency K
  best_aic <- Inf
  best_k <- NA
  best_beta <- NULL

  for (k in k_candidates) {

    # Sin and cos founction
    sin_term <- sin(2 * pi * k * t_index / N_total)
    cos_term <- cos(2 * pi * k * t_index / N_total)
    X_aug <- cbind(1, x_trim, sin_term, cos_term, dx_matrix)

    # OLS Estimation: (X'X)^(-1) X'Y
    beta_hat <- solve(t(X_aug) %*% X_aug) %*% t(X_aug) %*% y_trim

    # SSR calculation for AIC
    residuals <- y_trim - (X_aug %*% beta_hat)
    ssr <- sum(residuals^2)
    p <- ncol(X_aug)

    # AIC calculation
    current_aic <- N_trim * log(ssr / N_trim) + 2 * p

    # To chose the optimal k using AIC
    if (current_aic < best_aic) {
      best_aic <- current_aic
      best_k <- k

      # To name rows
      rownames(beta_hat)[1:4] <- c("Intercept", "X_Coef", "Gamma_1 (Sin)", "Gamma_2 (Cos)")
      best_beta <- beta_hat
    }
  }

# No4. To list the output
  return(list(
    Optimal_k = best_k,
    AIC = best_aic,
    Fourier_DOLS_Estimates = best_beta[1:4, 1, drop = FALSE]
  ))
}
