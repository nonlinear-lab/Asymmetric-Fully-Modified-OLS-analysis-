dols_manual <- function(y, x, n_leads, n_lags) {
  # Ensure inputs are matrices
  y <- as.matrix(y)
  x <- as.matrix(x)
  T <- nrow(x)

  # 1. First difference of x
  dx <- diff(x)

  # 2. Create leads and lags matrix for dx
  # We lose 1 observation due to differencing, plus leads and lags
  valid_idx <- (n_lags + 2):(T - n_leads)

  # Trim y and x to match the valid indices
  y_trim <- y[valid_idx, , drop = FALSE]
  x_trim <- x[valid_idx, , drop = FALSE]

  # Initialize the augmented regressor matrix with the base x variable
  X_aug <- x_trim

  # Populate leads and lags of dx
  for (j in (-n_leads):n_lags) {
    # If j is negative, it's a lead. If positive, it's a lag.
    dx_shifted <- dx[(valid_idx - 1 - j), , drop = FALSE]
    X_aug <- cbind(X_aug, dx_shifted)
  }

  # 3. Add an intercept
  X_aug <- cbind(1, X_aug)

  # 4. Standard OLS Estimation: (X'X)^(-1) X'Y
  beta_hat <- solve(t(X_aug) %*% X_aug) %*% t(X_aug) %*% y_trim

  # The first coefficient is the intercept, the second is our cointegrating vector
  return(beta_hat[2, 1]) 
}