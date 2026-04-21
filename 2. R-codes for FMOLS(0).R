#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: set a dorking directory and library
#++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/R/fmols")

# load the package
# install.packages("cointReg")
library(cointReg)
mout<- matrix(data=NA,nrow=4,ncol=1)
#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: load data
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
  #++++++++++++++++++++++++++++++++++++++++++++++++
  # Step 3: DOLS analysis
  #++++++++++++++++++++++++++++++++++++++++++++++++
  # 2. DOLS Estimation
  # Create a constant intercept
  intercept <- rep(1, length(y))

  # We use the cointRegFM function. You need to specify kernel and bandwith.
  # Here we use 2 leads and 2 lags as an example.
  fmols_model <- cointRegFM(y = y, x = x, deter = intercept, kernel= "ba", bandwith = 2)
  fmols   <- fmols_model$theta[2]
  mout[j,1]   <- fmols
  print(fmols_model)
}
print(mout)
