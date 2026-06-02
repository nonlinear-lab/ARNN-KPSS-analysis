# =====================================================================
# MONTE CARLO SIMULATION: SIZE & POWER ANALYSIS + GRAPH GENERATION
# Includes: KPSS, FKPSS, KPSS-SB, FKPSS-SB, ARNN-KPSS (p=1)
# =====================================================================

# Set your Working Directory (WD) where the graphs will be saved
setwd("C:/R/arnn")

library(ggplot2)
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
# PART 2: Test Statistics (5 Models)
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
# PART 3: Simulation Execution Wrapper
# ---------------------------------------------------------------------
run_full_simulation <- function(n_sim = 1000, n_obs = 100, lambda = 0, rho = 0, omega = 0.125, k = 1, cvs) {
  l0 <- 0
  l4 <- floor(4 * (n_obs / 100)^0.25)
  l12 <- floor(12 * (n_obs / 100)^0.25)
  
  rejections <- matrix(0, nrow = 5, ncol = 3)
  rownames(rejections) <- c("KPSS", "FKPSS", "KPSS_SB", "FKPSS_SB", "ARNN_KPSS")
  colnames(rejections) <- c("L0", "L4", "L12")
  
  for (i in 1:n_sim) {
    y <- simulate_kpss_dgp(n_obs, lambda, rho)
    
    if (kpss_test(y, l0) > cvs$kpss) rejections["KPSS", "L0"] <- rejections["KPSS", "L0"] + 1
    if (kpss_test(y, l4) > cvs$kpss) rejections["KPSS", "L4"] <- rejections["KPSS", "L4"] + 1
    if (kpss_test(y, l12) > cvs$kpss) rejections["KPSS", "L12"] <- rejections["KPSS", "L12"] + 1
    
    if (fkpss_test(y, k, l0) > cvs$fkpss) rejections["FKPSS", "L0"] <- rejections["FKPSS", "L0"] + 1
    if (fkpss_test(y, k, l4) > cvs$fkpss) rejections["FKPSS", "L4"] <- rejections["FKPSS", "L4"] + 1
    if (fkpss_test(y, k, l12) > cvs$fkpss) rejections["FKPSS", "L12"] <- rejections["FKPSS", "L12"] + 1
    
    if (kpss_sb_test(y, omega, l0) > cvs$kpss_sb) rejections["KPSS_SB", "L0"] <- rejections["KPSS_SB", "L0"] + 1
    if (kpss_sb_test(y, omega, l4) > cvs$kpss_sb) rejections["KPSS_SB", "L4"] <- rejections["KPSS_SB", "L4"] + 1
    if (kpss_sb_test(y, omega, l12) > cvs$kpss_sb) rejections["KPSS_SB", "L12"] <- rejections["KPSS_SB", "L12"] + 1
    
    if (fkpss_sb_test(y, k, omega, l0) > cvs$fkpss_sb) rejections["FKPSS_SB", "L0"] <- rejections["FKPSS_SB", "L0"] + 1
    if (fkpss_sb_test(y, k, omega, l4) > cvs$fkpss_sb) rejections["FKPSS_SB", "L4"] <- rejections["FKPSS_SB", "L4"] + 1
    if (fkpss_sb_test(y, k, omega, l12) > cvs$fkpss_sb) rejections["FKPSS_SB", "L12"] <- rejections["FKPSS_SB", "L12"] + 1
    
    if (arnn_kpss_test(y, l0) > cvs$arnn) rejections["ARNN_KPSS", "L0"] <- rejections["ARNN_KPSS", "L0"] + 1
    if (arnn_kpss_test(y, l4) > cvs$arnn) rejections["ARNN_KPSS", "L4"] <- rejections["ARNN_KPSS", "L4"] + 1
    if (arnn_kpss_test(y, l12) > cvs$arnn) rejections["ARNN_KPSS", "L12"] <- rejections["ARNN_KPSS", "L12"] + 1
  }
  
  return(rejections / n_sim)
}

# ---------------------------------------------------------------------
# PART 4: Execution Setup and ggplot2 Graph Generation
# ---------------------------------------------------------------------
set.seed(123)
n_values <- c(50, 100, 500, 1000)
bw_to_plot <- "L12" # Bandwidth to use for plotting

# Helper function to get correct 5% CVs for each sample size
get_critical_values <- function(n) {
  if (n == 50) return(list(kpss=0.143, fkpss=0.058, kpss_sb=0.116, fkpss_sb=0.048, arnn=0.141))
  if (n == 100) return(list(kpss=0.147, fkpss=0.053, kpss_sb=0.118, fkpss_sb=0.046, arnn=0.145))
  if (n == 500) return(list(kpss=0.145, fkpss=0.052, kpss_sb=0.115, fkpss_sb=0.046, arnn=0.145))
  if (n == 1000) return(list(kpss=0.146, fkpss=0.053, kpss_sb=0.114, fkpss_sb=0.046, arnn=0.144))
  return(list(kpss=0.146, fkpss=0.053, kpss_sb=0.115, fkpss_sb=0.046, arnn=0.145))
}

# Visual Mapping for ggplot2 
test_labels <- c("KPSS", "FKPSS", "KPSS_SB", "FKPSS_SB", "ARNN_KPSS")
line_types  <- c("dashed", "dotted", "longdash", "dotdash", "solid")
point_shapes<- c(15, 17, 4, 18, 16)

rho_values <- c(0, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99)
lambda_values <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100)

# Master Loop across all sample sizes
for (n_val in n_values) {
  
  cat(sprintf("\n\n========================================================\n"))
  cat(sprintf(" RUNNING ANALYSIS FOR SAMPLE SIZE N = %d \n", n_val))
  cat(sprintf("========================================================\n"))
  
  critical_values <- get_critical_values(n_val)
  
  # ---------------- A. SIZE ANALYSIS ----------------
  cat(sprintf("\n--- A. SIZE ANALYSIS (Lambda = 0) ---\n"))
  size_data <- data.frame()
  
  for (r in rho_values) {
    cat(sprintf("Size Results for Rho = %.2f ...\n", r))
    res_size <- run_full_simulation(n_sim = 1000, n_obs = n_val, lambda = 0, rho = r, 
                                    omega = 0.125, k = 1, cvs = critical_values)
    print(round(res_size, 3))
    
    res_extracted <- res_size[, bw_to_plot]
    df <- data.frame(Rho = r, Test = factor(names(res_extracted), levels = test_labels), RejectionRate = res_extracted)
    size_data <- rbind(size_data, df)
  }
  
  p_size <- ggplot(size_data, aes(x = Rho, y = RejectionRate, group = Test, 
                                  linetype = Test, shape = Test)) +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = rho_values) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, 0.2)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = point_shapes) +
    theme_bw(base_size = 14) +
    labs(title = sprintf("Empirical Size (Bandwidth %s, N = %d)", bw_to_plot, n_val),
         x = "Autocorrelation Parameter (\u03C1)",
         y = "Empirical Size") +
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # ---------------- B. POWER ANALYSIS ----------------
  cat(sprintf("\n--- B. POWER ANALYSIS (Rho = 0) ---\n"))
  power_data <- data.frame()
  
  for (lam in lambda_values) {
    cat(sprintf("Power Results for Lambda = %f ...\n", lam))
    res_power <- run_full_simulation(n_sim = 1000, n_obs = n_val, lambda = lam, rho = 0, 
                                     omega = 0.125, k = 1, cvs = critical_values)
    print(round(res_power, 3))
    
    res_extracted <- res_power[, bw_to_plot]
    df <- data.frame(Lambda = lam, Test = factor(names(res_extracted), levels = test_labels), RejectionRate = res_extracted)
    power_data <- rbind(power_data, df)
  }
  
  p_power <- ggplot(power_data, aes(x = Lambda, y = RejectionRate, group = Test, 
                                    linetype = Test, shape = Test)) +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_log10(breaks = lambda_values, 
                  labels = c("1e-4", "1e-3", "0.01", "0.1", "1", "10", "100")) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, 0.2)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = point_shapes) +
    theme_bw(base_size = 14) +
    labs(title = sprintf("Empirical Power (Bandwidth %s, N = %d)", bw_to_plot, n_val),
         x = expression("Variance Ratio: " * lambda * " = " * sigma[u]^2 / sigma[epsilon]^2),
         y = "Empirical Power") +
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Print the plots to the viewer
  print(p_size)
  print(p_power)
  
  # Automatically save the 8 plots directly to the Working Directory as high-quality TIFF files
  ggsave(sprintf("Size_%s_N%d.tiff", bw_to_plot, n_val), plot = p_size, width = 9, height = 6, dpi = 300, device = "tiff")
  ggsave(sprintf("Power_%s_N%d.tiff", bw_to_plot, n_val), plot = p_power, width = 9, height = 6, dpi = 300, device = "tiff")
}