#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: set a dorking directory and library
#++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/R/fmols")

# load the package
# install.packages("cointReg")
library(cointReg)
# mout<-  matrix(NA, nrow=4, ncol=1)
mout <- matrix(NA, nrow=4, ncol=2,
               dimnames = list(paste0("Model", 1:4), c("beta+", "beta-")))
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

  theta <- solve(XX) %*% (Yplus_X - correction)

  # Return the stochastic slope (beta)
  return(list(beta_pos = theta[2, 1], beta_neg = theta[3, 1]))  #
}
#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: load data and run
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(4).csv")
my_data <- as.matrix(my_data)
n       <- nrow(my_data)

for (j in 1:4) {
  if (j == 1) {
    y <- as.matrix(my_data[2:n, 3])
    x <- as.matrix(my_data[2:n, 2])
  } else if (j == 2) {
    y <- as.matrix(my_data[2:n, 4])
    x <- as.matrix(my_data[2:n, 2])
  } else if (j == 3) {
    set.seed(42)
    y <- as.matrix(my_data[2:n, 5])
    x <- as.matrix(my_data[2:n, 2])
  } else {
    set.seed(42)
    y <- as.matrix(my_data[2:n, 6])
    x <- as.matrix(my_data[2:n, 2])
  }
  # Call the function using your 'y' and 'x' data, specifying leads and lags
  fmols <-manual_fmols (y, x, bandwidth = 2)
  mout[j, 1] <- fmols$beta_pos
  mout[j, 2] <- fmols$beta_neg
  print(fmols)
}
print(mout)
