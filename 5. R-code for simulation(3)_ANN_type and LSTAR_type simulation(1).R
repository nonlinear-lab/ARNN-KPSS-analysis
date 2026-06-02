# =====================================================================
# MONTE CARLO SIMULATION SCRIPT: ARNN-KPSS POWER ANALYSIS
# Includes True ARNN and LSTAR Data Generating Processes
# =====================================================================

# ---------------------------------------------------------------------
# PART 1: Data Generating Processes (DGPs)
# ---------------------------------------------------------------------

# 1A. True Autoregressive Neural Network (ARNN) DGP
simulate_true_arnn <- function(n, lambda, gamma_0 = 0, gamma_1 = 1, 
                               theta = 0.5, rho = 0.5, sd_eps = 1) {
  var_eps <- sd_eps^2
  var_u <- lambda * var_eps
  sd_u <- sqrt(var_u)
  
  y <- numeric(n)
  eps <- rnorm(n, mean = 0, sd = sd_eps)
  u <- rnorm(n, mean = 0, sd = sd_u)
  r <- cumsum(u) 
  
  y[1] <- eps[1] + r[1]
  for (t in 2:n) {
    # Logistic activation function bypassing Taylor truncation
    activation <- 1 / (1 + exp(-(gamma_0 + gamma_1 * y[t-1])))
    y[t] <- rho * y[t-1] + theta * activation + r[t] + eps[t]
  }
  return(y)
}

# 1B. Logistic Smooth Transition Autoregressive (LSTAR) DGP
simulate_lstar <- function(n, lambda, alpha_1 = 0.6, alpha_2 = -0.4, gamma = 2, c = 0, sd_eps = 1) {
  var_eps <- sd_eps^2
  var_u <- lambda * var_eps
  sd_u <- sqrt(var_u)
  
  y <- numeric(n)
  eps <- rnorm(n, mean = 0, sd = sd_eps)
  u <- rnorm(n, mean = 0, sd = sd_u)
  r <- cumsum(u) 
  
  y[1] <- eps[1] + r[1]
  for (t in 2:n) {
    # Smooth transition mechanism
    transition <- 1 / (1 + exp(-gamma * (y[t-1] - c)))
    y[t] <- (alpha_1 * y[t-1]) + (alpha_2 * y[t-1]) * transition + r[t] + eps[t]
  }
  return(y)
}


# ---------------------------------------------------------------------
# PART 2: Helper Functions & Test Statistics
# ---------------------------------------------------------------------

# 2A. Long-Run Variance (Bartlett Window)
calc_lrv <- function(e, l) {
  n <- length(e)
  var_e <- sum(e^2) / n
  
  if (l > 0) {
    cov_e <- 0
    for (s in 1:l) {
      weight <- 1 - (s / (l + 1))
      cov_s <- sum(e[(s + 1):n] * e[1:(n - s)]) / n
      cov_e <- cov_e + weight * cov_s
    }
    lrv <- var_e + 2 * cov_e
  } else {
    lrv <- var_e
  }
  return(lrv)
}

# 2B. Standard KPSS Test
kpss_test <- function(y, l) {
  n <- length(y)
  t <- 1:n
  mod <- lm(y ~ t)
  e <- residuals(mod)
  S <- cumsum(e)
  eta <- sum(S^2) / (n^2 * calc_lrv(e, l))
  return(eta)
}

# 2C. Fourier KPSS (FKPSS) Test
fkpss_test <- function(y, k = 1, l) {
  n <- length(y)
  t <- 1:n
  sin_t <- sin(2 * pi * k * t / n)
  cos_t <- cos(2 * pi * k * t / n)
  mod <- lm(y ~ t + sin_t + cos_t)
  e <- residuals(mod)
  S <- cumsum(e)
  eta <- sum(S^2) / (n^2 * calc_lrv(e, l))
  return(eta)
}

# 2D. ARNN-KPSS Test
arnn_kpss_test <- function(y, l) {
  n <- length(y)
  y_curr <- y[2:n]
  t_trend <- 2:n
  y_lag1 <- y[1:(n - 1)]
  y_lag1_sq <- y_lag1^2
  y_lag1_cu <- y_lag1^3
  
  mod <- lm(y_curr ~ t_trend + y_lag1 + y_lag1_sq + y_lag1_cu)
  e <- residuals(mod)
  S <- cumsum(e)
  n_adj <- length(e)
  eta <- sum(S^2) / (n_adj^2 * calc_lrv(e, l))
  return(eta)
}


# ---------------------------------------------------------------------
# PART 3: Monte Carlo Simulation Engine
# ---------------------------------------------------------------------

run_power_simulation <- function(n_sim = 1000, n_obs = 100, lambda = 0.1, 
                                 cv_arnn, cv_kpss, cv_fkpss, dgp_type = "ARNN") {
  # Calculate bandwidths
  l0 <- 0
  l4 <- floor(4 * (n_obs / 100)^0.25)
  l12 <- floor(12 * (n_obs / 100)^0.25)
  
  # Initialize rejection counters
  rej <- list(arnn = c(l0=0, l4=0, l12=0), 
              kpss = c(l0=0, l4=0, l12=0), 
              fkpss = c(l0=0, l4=0, l12=0))
  
  for (i in 1:n_sim) {
    # Select DGP
    if (dgp_type == "ARNN") {
      y <- simulate_true_arnn(n_obs, lambda)
    } else {
      y <- simulate_lstar(n_obs, lambda)
    }
    
    # ARNN-KPSS
    if(arnn_kpss_test(y, l0) > cv_arnn) rej$arnn['l0'] <- rej$arnn['l0'] + 1
    if(arnn_kpss_test(y, l4) > cv_arnn) rej$arnn['l4'] <- rej$arnn['l4'] + 1
    if(arnn_kpss_test(y, l12) > cv_arnn) rej$arnn['l12'] <- rej$arnn['l12'] + 1
    
    # KPSS
    if(kpss_test(y, l0) > cv_kpss) rej$kpss['l0'] <- rej$kpss['l0'] + 1
    if(kpss_test(y, l4) > cv_kpss) rej$kpss['l4'] <- rej$kpss['l4'] + 1
    if(kpss_test(y, l12) > cv_kpss) rej$kpss['l12'] <- rej$kpss['l12'] + 1
    
    # FKPSS (k=1)
    if(fkpss_test(y, 1, l0) > cv_fkpss) rej$fkpss['l0'] <- rej$fkpss['l0'] + 1
    if(fkpss_test(y, 1, l4) > cv_fkpss) rej$fkpss['l4'] <- rej$fkpss['l4'] + 1
    if(fkpss_test(y, 1, l12) > cv_fkpss) rej$fkpss['l12'] <- rej$fkpss['l12'] + 1
  }
  
  # Return power as a formatted row
  return(data.frame(
    Lambda = lambda,
    ARNN_L0 = rej$arnn['l0']/n_sim, ARNN_L4 = rej$arnn['l4']/n_sim, ARNN_L12 = rej$arnn['l12']/n_sim,
    KPSS_L0 = rej$kpss['l0']/n_sim, KPSS_L4 = rej$kpss['l4']/n_sim, KPSS_L12 = rej$kpss['l12']/n_sim,
    FKPSS_L0 = rej$fkpss['l0']/n_sim, FKPSS_L4 = rej$fkpss['l4']/n_sim, FKPSS_L12 = rej$fkpss['l12']/n_sim,
    row.names = NULL
  ))
}


# ---------------------------------------------------------------------
# PART 4: Execution Block (Generating the Tables)
# ---------------------------------------------------------------------

set.seed(123)

# 1. Setup Parameters (Example for n = 100)
n_obs <- 100
n_simulations <- 1000 
lambda_values <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100) 

# 2. Extract 5% Critical Values for n=100
cv_arnn_n100 <- 0.062 
cv_kpss_n100 <- 0.155  
cv_fkpss_n100 <- 0.056

cat("========================================================\n")
cat("RUNNING SIMULATIONS FOR TRUE ARNN DGP (n = 100)\n")
cat("========================================================\n")

# Run loop over lambdas for True ARNN
results_arnn_dgp <- do.call(rbind, lapply(lambda_values, function(lam) {
  run_power_simulation(n_sim = n_simulations, n_obs = n_obs, lambda = lam, 
                       cv_arnn = cv_arnn_n100, cv_kpss = cv_kpss_n100, 
                       cv_fkpss = cv_fkpss_n100, dgp_type = "ARNN")
}))
print(results_arnn_dgp)

cat("\n========================================================\n")
cat("RUNNING SIMULATIONS FOR LSTAR DGP (n = 100)\n")
cat("========================================================\n")

# Run loop over lambdas for LSTAR
results_lstar_dgp <- do.call(rbind, lapply(lambda_values, function(lam) {
  run_power_simulation(n_sim = n_simulations, n_obs = n_obs, lambda = lam, 
                       cv_arnn = cv_arnn_n100, cv_kpss = cv_kpss_n100, 
                       cv_fkpss = cv_fkpss_n100, dgp_type = "LSTAR")
}))
print(results_lstar_dgp)