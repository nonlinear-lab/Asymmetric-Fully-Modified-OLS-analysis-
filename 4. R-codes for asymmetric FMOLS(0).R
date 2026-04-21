#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: set a working directory and library
#++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/R/fmols")
library(cointReg)

mout <- matrix(data=NA, nrow=4, ncol=2)
#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: load data
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(4).csv")
my_data <- as.matrix(my_data)
n       <- nrow(my_data)

for (j in 1:4) {
  # Standardize indexing to 2:T (199 observations)
  if (j == 1) {
    y <- as.matrix(my_data[2:n, 3])
    x <- as.matrix(my_data[2:n, 2])
    k_val <- 4
  } else if (j == 2) {
    y <- as.matrix(my_data[2:n, 4])
    x <- as.matrix(my_data[2:n, 2])
    k_val <- 2
  } else if (j == 3) {
    y <- as.matrix(my_data[2:n, 5])
    x <- as.matrix(my_data[2:n, 2])
    k_val <- 1
  } else {
    y <- as.matrix(my_data[2:n, 6])
    x <- as.matrix(my_data[2:n, 2])
    k_val <- 1
  }

  #++++++++++++++++++++++++++++++++++++++++++++++++
  # Step 3: FMOLS analysis
  #++++++++++++++++++++++++++++++++++++++++++++++++
  T_total <- nrow(y)
  t_index <- 1:T_total
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
  y_trim     <- y[-1, , drop = FALSE]
  x_asym     <- cbind(x_pos, x_neg)
  deter_mat  <- matrix(1, nrow(y_trim), 1)

    # FMOLS Estimation
  fmols_model <- cointRegFM(y = y_trim,
                            x = x_asym,
                            deter = deter_mat,
                            kernel = "ba",
                            bandwidth = 2)

  # FIX 3: theta[4] is the X coefficient
  # (theta[1]=Int, [2]=Sin, [3]=Cos, [4]=X)
  mout[j,1]  <- fmols_model$theta[2]
  mout[j,2]  <- fmols_model$theta[3]

  cat("\n--- Scenario", j, " (k =", k_val, ") ---\n")
  print(fmols_model)
}

cat("\nFinal Results (X-Slopes):\n")
print(mout)
