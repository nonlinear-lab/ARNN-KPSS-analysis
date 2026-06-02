# =====================================================================
# MONTE CARLO SIMULATION: EMPIRICAL DISTRIBUTION OF KPSS-TYPE TESTS
# Generates Tables and ECDF Plots for N = 50, 100, 500, 1000
# =====================================================================

# Install necessary packages if you don't have them:
# install.packages(c("ggplot2", "tidyr", "dplyr"))

library(ggplot2)
library(tidyr)
library(dplyr)

# Optional: Set your Working Directory (WD) where the graphs will be saved
setwd("C:/R/arnn") 

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
# Instead of returning just the matrix, we return all 1000 raw statistics 
# so we can build both the table and the ECDF graph.
simulate_empirical_stats <- function(n_sim = 1000, n_obs = 100, k = 1, omega = 0.125) {
  
  stats <- list(
    KPSS      = numeric(n_sim),
    FKPSS     = numeric(n_sim),
    KPSS_SB   = numeric(n_sim),
    FKPSS_SB  = numeric(n_sim),
    ARNN_KPSS = numeric(n_sim)
  )
  l0 <- 0 
  
  for (i in 1:n_sim) {
    y <- simulate_kpss_dgp(n_obs, lambda = 0, rho = 0)
    
    stats$KPSS[i]      <- kpss_test(y, l0)
    stats$FKPSS[i]     <- fkpss_test(y, k, l0)
    stats$KPSS_SB[i]   <- kpss_sb_test(y, omega, l0)
    stats$FKPSS_SB[i]  <- fkpss_sb_test(y, k, omega, l0)
    stats$ARNN_KPSS[i] <- arnn_kpss_test(y, l0)
  }
  return(stats)
}

# ---------------------------------------------------------------------
# PART 5: Execution Setup and Graph Generation
# ---------------------------------------------------------------------
set.seed(123)

n_values <- c(50, 100, 500, 1000)
k_target <- 1
omega_target <- 0.125
n_simulations <- 1000

# Plotting parameters mapping to your image style
plot_probs <- c(0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99)
test_labels <- c("ARNN-KPSS", "KPSS", "FKPSS (k=1)", "KPSS-SB", "FKPSS-SB (k=1, w=0.125)")
line_types  <- c("solid", "dashed", "dotted", "longdash", "dotdash")
point_shapes<- c(16, 15, 17, 4, 18) # 16=Circle, 15=Square, 17=Triangle, 4=Cross, 18=Diamond

cat("========================================================\n")
cat(" ESTIMATING EMPIRICAL CRITICAL VALUES & PLOTTING \n")
cat("========================================================\n")

for (n in n_values) {
  cat(sprintf("\n--- Simulating for n = %d, k = %d, omega = %.3f ---\n", 
              n, k_target, omega_target))
  
  # 1. Run Simulation
  stats <- simulate_empirical_stats(n_sim = n_simulations, n_obs = n, 
                                    k = k_target, omega = omega_target)
  
  # 2. Print the Upper-Tail Critical Value Matrix to Console
  cv_percentiles <- c(0.90, 0.95, 0.99)
  cv_matrix <- rbind(
    KPSS      = quantile(stats$KPSS, cv_percentiles),
    FKPSS     = quantile(stats$FKPSS, cv_percentiles),
    KPSS_SB   = quantile(stats$KPSS_SB, cv_percentiles),
    FKPSS_SB  = quantile(stats$FKPSS_SB, cv_percentiles),
    ARNN_KPSS = quantile(stats$ARNN_KPSS, cv_percentiles)
  )
  colnames(cv_matrix) <- c("10%", "5%", "1%")
  print(round(cv_matrix, 3))
  
  # 3. Prepare Data for ECDF Graph
  plot_df <- data.frame(
    CumulativeProbability = rep(plot_probs, 5),
    CriticalValue = c(
      quantile(stats$ARNN_KPSS, plot_probs),
      quantile(stats$KPSS, plot_probs),
      quantile(stats$FKPSS, plot_probs),
      quantile(stats$KPSS_SB, plot_probs),
      quantile(stats$FKPSS_SB, plot_probs)
    ),
    Test = factor(rep(test_labels, each = length(plot_probs)), levels = test_labels)
  )
  
  # 4. Generate ggplot
  p <- ggplot(plot_df, aes(x = CriticalValue, y = CumulativeProbability, 
                           group = Test, linetype = Test, shape = Test)) +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(limits = c(0.00, 0.26), breaks = seq(0.00, 0.25, 0.05)) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, 0.2)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = point_shapes) +
    theme_bw(base_size = 14) +
    labs(title = sprintf("Empirical Cumulative Distribution Function (N = %d)", n),
         x = "Critical Value",
         y = "Cumulative Probability") +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 5. Print to Viewer and Save
  print(p)
  filename <- sprintf("Dist_ECDF_N%d.tiff", n)
  ggsave(filename, plot = p, width = 9, height = 6, dpi = 300, device = "tiff")
  
  cat(sprintf("-> Saved graph as '%s' in %s\n", filename, getwd()))
}