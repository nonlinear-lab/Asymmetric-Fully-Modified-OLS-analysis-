#+++++++++++++++++++++++++++++++++
# Step 0: set a dorking directory
#+++++++++++++++++++++++++++++++
setwd("C:/R/fmols")

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: to create a loop and define y and x
#++++++++++++++++++++++++++++++++++++++++++++++
T <- 200
set.seed(42)
u_x <- rnorm(T, mean = 0, sd = 1)
x <- cumsum(u_x)

# 1. Generate y as four version
for (j in 1:6) {
  if (j == 1) {
    set.seed(99)
    u_y <- rnorm(T, mean = 0, sd = 1)
    y1 <- 10 + 2.5 * x + u_y
  } else if (j == 2) {
    set.seed(99)
    u_y <- arima.sim(model = list(ar = 0), n = T, sd = 1)
    y2 <- 10 + 2.5 * x + u_y
  } else if (j == 3) {
    set.seed(99)
    u_y <- arima.sim(model = list(ar = 0.8), n = T, sd = 1)
    y3 <- 10 + 2.5 * x + u_y
  } else if (j == 4) {
    set.seed(99)
    u_y <- arima.sim(model = list(ar = 0.95), n = T, sd = 1)
    y4 <- 10 + 2.5 * x + u_y
  } else if (j == 5) {
    set.seed(99)
    u_y <- arima.sim(model = list(order = c(0,1,0)), n = T, sd = 1)
    y5 <- 10 + 2.5 * x + u_y
  } else {
    set.seed(99)
    u_y <- cumsum(rnorm(T, mean = 0, sd = 1))
    y6 <- 10 + 2.5 * x + u_y
  }
}

#+++++++++++++++++++++++++++++++++
# Step 2:  Generate data & Adjust Margins
#+++++++++++++++++++++++++++++++
y_limits <- range(c(x, y1, y2, y3, y4))

# 1. Set a wide right margin to make room (10 lines of space)
par(mar = c(5, 4, 4, 10))

# 2. Plot the data
plot(x, type="l", ylim=y_limits, col="black", ylab="Values", xlab="Time", main="Simulation")
lines(y1, col="blue")
lines(y2, col="yellow")
lines(y3, col="green")
lines(y4, col="purple")
lines(y5, col="pink")
lines(y6, col="red")

#+++++++++++++++++++++++++++++++++
# Step 3: Create legend OUTSIDE the graph
#+++++++++++++++++++++++++++++++
# 3. Use exact coordinates: Start drawing at X=205, Y=Maximum height
legend(x = 205, y = max(y_limits), xpd = TRUE,
       legend=c("X (Random Walk)", "Strong (white noise)", "Strong (AR=0)", "Weak (AR=0.8)", "Very Weak (AR=0.95)", "None (AR=1.0)", "None (random walk)"),
       col=c("black", "blue", "yellow", "green", "purple", "pink","red"),
       lty=1, lwd=2, cex=0.55, bty="n")

# Reset the margins back to normal when finished
par(mar = c(5, 4, 4, 2))

#+++++++++++++++++++++++++++++++++
# Step 4: Save the Data
#+++++++++++++++++++++++++++++++
# Combine x and all four y variables into a single data frame
master_data <- data.frame(
  Time_Index = 1:T,
  X_Random_Walk = x,
  Y1_Strong_Coint = y1,
  Y2_Weak_Coint = y2,
  Y3_Very_Weak_Coint = y3,
  Y4_No_Coint = y4
)
# Save the combined data to your hard drive
write.csv(master_data, file = "data(5).csv", row.names = FALSE)
