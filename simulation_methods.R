#simulation_analysis <- function(scenario, T_val=200, reps=1000, true_beta=2.5,
                                # bandwidth=8, n_leads=2, n_lags= 2) {
simulation_analysis <- function(scenario,T_val,reps,true_beta,bandwidth,n_leads,n_lags) {

  results_slope <- matrix(NA, nrow = reps, ncol = 9)
  colnames(results_slope) <- c(
    "OLS", "DOLS", "Fourier_DOLS", "FMOLS", "Fourier_FMOLS",
    "Asym_DOLS_plus", "Asym_DOLS_minus",
    "Asym_FMOLS_plus", "Asym_FMOLS_minus"
  )

  for (i in 1:reps) {
    x <- cumsum(rnorm(T_val))
    u <- switch(scenario,
      "1" = rnorm(T_val),
      "2" = arima.sim(list(ar = 0.80), T_val),
      "3" = arima.sim(list(ar = 0.95), T_val),
      "4" = cumsum(rnorm(T_val))
    )
    y <- 10 + true_beta * x + u

    results_slope[i, 1] <- coef(lm(y ~ x))[2]
    results_slope[i, 2] <- dols_manual(y, x, n_leads, n_lags)
    results_slope[i, 3] <- fourier_dols_manual(y, x, n_leads, n_lags)$Fourier_DOLS_Estimates["X_Coef", 1]
    results_slope[i, 4] <- fmols_manual(y, x, bandwidth)
    results_slope[i, 5] <- fourier_fmols_manual(y, x, bandwidth)$Fourier_FMOLS_Estimates["X_Coef", 1]
    asym_res1 <- asymmetric_dols_manual(y, x, n_leads, n_lags)
    results_slope[i, 6] <- asym_res1$beta_plus
    results_slope[i, 7] <- asym_res1$beta_minus
    asym_res2 <- asymmetric_fmols_manual(y, x, bandwidth)
    results_slope[i, 8] <- asym_res2$beta_pos
    results_slope[i, 9] <- asym_res2$beta_neg
    if (i %% 100 == 0) cat("Repetition", i, "Complete\n")
  }
  summary_stats <- t(apply(results_slope, 2, function(col) {  # ✅ col used correctly
    c(Mean = mean(col),
      Bias = mean(col) - true_beta,
      SD   = sd(col),
      RMSE = sqrt(mean((col - true_beta)^2)))
  }))
  return(summary_stats)
}
