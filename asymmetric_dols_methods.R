asymmetric_dols_manual <- function(y, x, n_leads, n_lags) {
  # Ensure inputs are matrices
  y <- as.matrix(y)
  x <- as.matrix(x)
  T <- nrow(x)

  # 1. Asymmetric Decomposition (Shin 2014)
  dx <- diff(x)
  dx_pos <- ifelse(dx > 0, dx, 0)
  dx_neg <- ifelse(dx < 0, dx, 0)

  # Calculate partial sums (x+ and x-)
  x_pos <- as.matrix(cumsum(c(0, dx_pos)))
  x_neg <- as.matrix(cumsum(c(0, dx_neg)))

  # 2. Alignment and Trimming
  # Valid indices accounting for differencing and leads/lags
  valid_idx <- (n_lags + 2):(T - n_leads)

  y_trim <- y[valid_idx, , drop = FALSE]
  # x_trim <- x[valid_idx, , drop = FALSE]
  x_pos_trim <- x_pos[valid_idx, , drop = FALSE]
  x_neg_trim <- x_neg[valid_idx, , drop = FALSE]

  # Initialize augmented matrix with the two long-run components
  # X_aug <- cbind(x_trim, x_pos_trim, x_neg_trim)
  X_aug <- cbind(x_pos_trim, x_neg_trim)

  # 3. Add Leads and Lags for BOTH positive and negative changes
  for (j in (-n_leads):n_lags) {
    # Lead/Lag for positive differences
    dx_pos_shifted <- dx_pos[(valid_idx - 1 - j), , drop = FALSE]
    # Lead/Lag for negative differences
    dx_neg_shifted <- dx_neg[(valid_idx - 1 - j), , drop = FALSE]

    X_aug <- cbind(X_aug, dx_pos_shifted, dx_neg_shifted)
  }

  # 4. Add Intercept
  X_aug <- cbind(1, X_aug)
  XtX_inv <- tryCatch(solve(t(X_aug) %*% X_aug), error = function(e) NULL)

  if (is.null(XtX_inv)) {
    return(list(beta_plus = NA, beta_minus = NA, se_plus = NA, se_minus = NA))
  }
  # 5. OLS Estimation: (X'X)^(-1) X'Y
  beta_hat <- solve(t(X_aug) %*% X_aug) %*% t(X_aug) %*% y_trim

  # Calculate Residuals and Conventional Variance
  u_hat <- y_trim - X_aug %*% beta_hat
  df <- nrow(X_aug) - ncol(X_aug)
  sigma_sq <- sum(u_hat^2) / df

  cov_matrix <- sigma_sq * XtX_inv
  se_hat <- sqrt(diag(cov_matrix))

  # Results: 1st=Intercept, 2nd=Beta+, 3rd=Beta-
  results <- list(
    beta_plus  = beta_hat[2, 1],
    beta_minus = beta_hat[3, 1],
    se_plus    = se_hat[2],
    se_minus   = se_hat[3]
  )
  return(results)
}
