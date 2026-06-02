# =====================================================================
# MONTE CARLO SIMULATION: EMPIRICAL DISTRIBUTION OF KPSS-TYPE TESTS
# Generates Full Table 1 for N = 50, 100, 500, 1000
# Includes combinations for k = 1, 2 and omega = 0.125, 0.375, 0.625, 0.875
# =====================================================================

library(tidyr)
library(dplyr)

# ---------------------------------------------------------------------
# PART 1: Data Generating Process (DGP)
# ---------------------------------------------------------------------
simulate_kpss_dgp <- function(n, lambda, rho = 0, mu = 0, beta = 0) {
  v <- rnorm(n, mean = 0, sd = 1)
  eps <- numeric(n)
  eps[1] <- v[1]
  if (n > 1) {
    for (t in 2:n) eps[t] <- rho * eps[t-1] + v[t]
  }
  var_eps <- 1 / (1 - rho^2) 
  var_u <- lambda * var_eps
  sd_u <- sqrt(var_u)
  u <- rnorm(n, mean = 0, sd = sd_u)
  r <- cumsum(u)
  t_seq <- 1:n
  return(mu + beta * t_seq + r + eps)
}

# ---------------------------------------------------------------------
# PART 2: Long-Run Variance (Bartlett Kernel)
# ---------------------------------------------------------------------
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

# ---------------------------------------------------------------------
# PART 3: Test Statistics
# ---------------------------------------------------------------------
kpss_test <- function(y, l) {
  t_seq <- 1:length(y)
  mod <- lm(y ~ t_seq)
  e <- residuals(mod)
  return(sum(cumsum(e)^2) / (length(e)^2 * calc_lrv(e, l)))
}

fkpss_test <- function(y, k, l) {
  n <- length(y)
  t_seq <- 1:n
  sin_t <- sin(2 * pi * k * t_seq / n)
  cos_t <- cos(2 * pi * k * t_seq / n)
  mod <- lm(y ~ t_seq + sin_t + cos_t)
  e <- residuals(mod)
  return(sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l)))
}

kpss_sb_test <- function(y, omega, l) {
  n <- length(y)
  t_seq <- 1:n
  n_B <- floor(omega * n)
  DU <- ifelse(t_seq > n_B, 1, 0)
  DT_B <- ifelse(t_seq == n_B, 1, 0)
  mod <- lm(y ~ t_seq + DU + DT_B)
  e <- residuals(mod)
  return(sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l)))
}

fkpss_sb_test <- function(y, k, omega, l) {
  n <- length(y)
  t_seq <- 1:n
  n_B <- floor(omega * n)
  sin_t <- sin(2 * pi * k * t_seq / n)
  cos_t <- cos(2 * pi * k * t_seq / n)
  DU <- ifelse(t_seq > n_B, 1, 0)
  DT_B <- ifelse(t_seq == n_B, 1, 0)
  mod <- lm(y ~ t_seq + sin_t + cos_t + DU + DT_B)
  e <- residuals(mod)
  return(sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l)))
}

arnn_kpss_test <- function(y, l) {
  n <- length(y)
  y_curr <- y[2:n]
  t_trend <- 2:n
  y_lag1 <- y[1:(n - 1)]
  mod <- lm(y_curr ~ t_trend + y_lag1 + I(y_lag1^2) + I(y_lag1^3))
  e <- residuals(mod)
  return(sum(cumsum(e)^2) / (length(e)^2 * calc_lrv(e, l)))
}

# ---------------------------------------------------------------------
# PART 4: Empirical Distribution Simulator
# ---------------------------------------------------------------------
simulate_empirical_stats <- function(n_sim = 1000, n_obs = 100) {
  
  # Initialize lists for all 16 requested test variations
  stats <- list(
    `ARNN-KPSS`               = numeric(n_sim),
    `KPSS`                    = numeric(n_sim),
    `FKPSS (k=1)`             = numeric(n_sim),
    `FKPSS (k=2)`             = numeric(n_sim),
    `KPSS-SB (w=0.125)`       = numeric(n_sim),
    `KPSS-SB (w=0.375)`       = numeric(n_sim),
    `KPSS-SB (w=0.625)`       = numeric(n_sim),
    `KPSS-SB (w=0.875)`       = numeric(n_sim),
    `FKPSS-SB (k=1, w=0.125)` = numeric(n_sim),
    `FKPSS-SB (k=1, w=0.375)` = numeric(n_sim),
    `FKPSS-SB (k=1, w=0.625)` = numeric(n_sim),
    `FKPSS-SB (k=1, w=0.875)` = numeric(n_sim),
    `FKPSS-SB (k=2, w=0.125)` = numeric(n_sim),
    `FKPSS-SB (k=2, w=0.375)` = numeric(n_sim),
    `FKPSS-SB (k=2, w=0.625)` = numeric(n_sim),
    `FKPSS-SB (k=2, w=0.875)` = numeric(n_sim)
  )
  l0 <- 0 
  
  for (i in 1:n_sim) {
    # Generate DGP under null (stationarity)
    y <- simulate_kpss_dgp(n_obs, lambda = 0, rho = 0)
    
    # 1. ARNN-KPSS & KPSS
    stats$`ARNN-KPSS`[i]         <- arnn_kpss_test(y, l0)
    stats$`KPSS`[i]              <- kpss_test(y, l0)
    
    # 2. FKPSS (k=1, 2)
    stats$`FKPSS (k=1)`[i]       <- fkpss_test(y, 1, l0)
    stats$`FKPSS (k=2)`[i]       <- fkpss_test(y, 2, l0)
    
    # 3. KPSS-SB (omega=0.125, 0.375, 0.625, 0.875)
    stats$`KPSS-SB (w=0.125)`[i] <- kpss_sb_test(y, 0.125, l0)
    stats$`KPSS-SB (w=0.375)`[i] <- kpss_sb_test(y, 0.375, l0)
    stats$`KPSS-SB (w=0.625)`[i] <- kpss_sb_test(y, 0.625, l0)
    stats$`KPSS-SB (w=0.875)`[i] <- kpss_sb_test(y, 0.875, l0)
    
    # 4. FKPSS-SB (k=1)
    stats$`FKPSS-SB (k=1, w=0.125)`[i] <- fkpss_sb_test(y, 1, 0.125, l0)
    stats$`FKPSS-SB (k=1, w=0.375)`[i] <- fkpss_sb_test(y, 1, 0.375, l0)
    stats$`FKPSS-SB (k=1, w=0.625)`[i] <- fkpss_sb_test(y, 1, 0.625, l0)
    stats$`FKPSS-SB (k=1, w=0.875)`[i] <- fkpss_sb_test(y, 1, 0.875, l0)
    
    # 5. FKPSS-SB (k=2)
    stats$`FKPSS-SB (k=2, w=0.125)`[i] <- fkpss_sb_test(y, 2, 0.125, l0)
    stats$`FKPSS-SB (k=2, w=0.375)`[i] <- fkpss_sb_test(y, 2, 0.375, l0)
    stats$`FKPSS-SB (k=2, w=0.625)`[i] <- fkpss_sb_test(y, 2, 0.625, l0)
    stats$`FKPSS-SB (k=2, w=0.875)`[i] <- fkpss_sb_test(y, 2, 0.875, l0)
  }
  return(stats)
}

# ---------------------------------------------------------------------
# PART 5: Execution Setup and Table Generation
# ---------------------------------------------------------------------
set.seed(123)

n_values <- c(50, 100, 500, 1000)
n_simulations <- 1000  # Note: Use 10000+ for publishable/robust results

cv_percentiles <- c(0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99)
col_names <- c("1%", "5%", "10%", "50%", "90%", "95%", "99%")

cat("========================================================\n")
cat(" TABLE 1: EMPIRICAL DISTRIBUTION OF ARNN-KPSS & KPSS TESTS \n")
cat("========================================================\n")

for (n in n_values) {
  cat(sprintf("\n--- Number of observations T = %d ---\n", n))
  
  # Run the full suite of simulations for the current n
  stats <- simulate_empirical_stats(n_sim = n_simulations, n_obs = n)
  
  # Calculate quantiles for every test in the list
  cv_matrix <- do.call(rbind, lapply(stats, function(x) quantile(x, cv_percentiles)))
  
  # Format table output to 3 decimal places
  colnames(cv_matrix) <- col_names
  print(round(cv_matrix, 3))
}