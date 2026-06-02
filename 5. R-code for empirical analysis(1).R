# =====================================================================
# EMPIRICAL ANALYSIS: KPSS, FKPSS, KPSS-SB, FKPSS-SB, ARNN-KPSS
# Automatically searches for optimal break dates (\omega) and frequencies (k)
# =====================================================================

library(dplyr)

# Optional: Set your Working Directory (WD) where the graphs will be saved
setwd("C:/R/arnn") 


# ---------------------------------------------------------------------
# 1. Long-Run Variance (Bartlett Kernel)
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
# 2. Test Statistics with Optimal Parameter Search (Minimizing SSR)
# ---------------------------------------------------------------------

# A. Standard KPSS
empirical_kpss <- function(y, l_vals) {
  t_seq <- 1:length(y)
  mod <- lm(y ~ t_seq)
  e <- residuals(mod)
  
  stats <- sapply(l_vals, function(l) sum(cumsum(e)^2) / (length(e)^2 * calc_lrv(e, l)))
  return(list(stats = stats))
}

# B. Fourier KPSS (Searches k = 1 and k = 2)
empirical_fkpss <- function(y, l_vals) {
  n <- length(y)
  t_seq <- 1:n
  
  best_ssr <- Inf
  best_mod <- NULL
  best_k <- NA
  
  for (k in 1:2) {
    sin_t <- sin(2 * pi * k * t_seq / n)
    cos_t <- cos(2 * pi * k * t_seq / n)
    mod <- lm(y ~ t_seq + sin_t + cos_t)
    ssr <- sum(residuals(mod)^2)
    
    if (ssr < best_ssr) {
      best_ssr <- ssr
      best_mod <- mod
      best_k <- k
    }
  }
  
  e <- residuals(best_mod)
  stats <- sapply(l_vals, function(l) sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l)))
  return(list(stats = stats, opt_k = best_k))
}

# C. KPSS with Structural Break (Searches omega between 0.15 and 0.85)
empirical_kpss_sb <- function(y, l_vals) {
  n <- length(y)
  t_seq <- 1:n
  
  trim <- 0.15
  lower_bound <- floor(trim * n)
  upper_bound <- floor((1 - trim) * n)
  
  best_ssr <- Inf
  best_mod <- NULL
  best_omega <- NA
  
  for (n_B in lower_bound:upper_bound) {
    DU <- ifelse(t_seq > n_B, 1, 0)
    DT_B <- ifelse(t_seq == n_B, 1, 0)
    mod <- lm(y ~ t_seq + DU + DT_B)
    ssr <- sum(residuals(mod)^2)
    
    if (ssr < best_ssr) {
      best_ssr <- ssr
      best_mod <- mod
      best_omega <- n_B / n
    }
  }
  
  e <- residuals(best_mod)
  stats <- sapply(l_vals, function(l) sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l)))
  return(list(stats = stats, opt_omega = round(best_omega, 2)))
}

# D. Fourier KPSS with Structural Break (Searches both k and omega)
empirical_fkpss_sb <- function(y, l_vals) {
  n <- length(y)
  t_seq <- 1:n
  
  trim <- 0.15
  lower_bound <- floor(trim * n)
  upper_bound <- floor((1 - trim) * n)
  
  best_ssr <- Inf
  best_mod <- NULL
  best_k <- NA
  best_omega <- NA
  
  for (k in 1:2) {
    sin_t <- sin(2 * pi * k * t_seq / n)
    cos_t <- cos(2 * pi * k * t_seq / n)
    
    for (n_B in lower_bound:upper_bound) {
      DU <- ifelse(t_seq > n_B, 1, 0)
      DT_B <- ifelse(t_seq == n_B, 1, 0)
      mod <- lm(y ~ t_seq + sin_t + cos_t + DU + DT_B)
      ssr <- sum(residuals(mod)^2)
      
      if (ssr < best_ssr) {
        best_ssr <- ssr
        best_mod <- mod
        best_k <- k
        best_omega <- n_B / n
      }
    }
  }
  
  e <- residuals(best_mod)
  stats <- sapply(l_vals, function(l) sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l)))
  return(list(stats = stats, opt_k = best_k, opt_omega = round(best_omega, 2)))
}

# E. ARNN-KPSS (p=1)
empirical_arnn <- function(y, l_vals) {
  n <- length(y)
  y_curr <- y[2:n]
  t_trend <- 2:n
  y_lag1 <- y[1:(n - 1)]
  
  mod <- lm(y_curr ~ t_trend + y_lag1 + I(y_lag1^2) + I(y_lag1^3))
  e <- residuals(mod)
  n_adj <- length(e)
  
  stats <- sapply(l_vals, function(l) sum(cumsum(e)^2) / (n_adj^2 * calc_lrv(e, l)))
  return(list(stats = stats))
}

# ---------------------------------------------------------------------
# 3. Main Execution Function
# ---------------------------------------------------------------------
run_country_analysis <- function(data_matrix) {
  countries <- colnames(data_matrix)
  n <- nrow(data_matrix)
  
  # Calculate Bandwidths
  l0 <- 0
  l4 <- floor(4 * (n / 100)^0.25)
  l12 <- floor(12 * (n / 100)^0.25)
  l_vals <- c(l0, l4, l12)
  
  results <- list()
  
  for (country in countries) {
    y <- na.omit(data_matrix[[country]]) # Drop NAs if any
    
    # Run all 5 tests
    res_kpss     <- empirical_kpss(y, l_vals)
    res_fkpss    <- empirical_fkpss(y, l_vals)
    res_kpss_sb  <- empirical_kpss_sb(y, l_vals)
    res_fkpss_sb <- empirical_fkpss_sb(y, l_vals)
    res_arnn     <- empirical_arnn(y, l_vals)
    
    # Format strings to match Tables 5 and 6 (e.g., "1.489 (1)", "0.922 [0.36]")
    fmt_kpss <- sprintf("%.3f", res_kpss$stats)
    fmt_fkpss <- sprintf("%.3f (%d)", res_fkpss$stats, res_fkpss$opt_k)
    fmt_kpss_sb <- sprintf("%.3f [%.2f]", res_kpss_sb$stats, res_kpss_sb$opt_omega)
    fmt_fkpss_sb <- sprintf("%.3f (%d)[%.2f]", res_fkpss_sb$stats, res_fkpss_sb$opt_k, res_fkpss_sb$opt_omega)
    fmt_arnn <- sprintf("%.3f", res_arnn$stats)
    
    country_df <- data.frame(
      Country = country,
      Bandwidth = c("L0", "L4", "L12"),
      KPSS = fmt_kpss,
      FKPSS = fmt_fkpss,
      KPSS_SB = fmt_kpss_sb,
      FKPSS_SB = fmt_fkpss_sb,
      ARNN_KPSS = fmt_arnn
    )
    
    results[[country]] <- country_df
  }
  
  final_table <- bind_rows(results)
  return(final_table)
}

# =====================================================================
# 4. LOAD YOUR DATA AND RUN (Example Placeholder)
# =====================================================================

# Instructions: Replace the lines below with your actual dataset.
# Your data should be a data frame where each column is a country's unemployment rate.

# EXAMPLE: 
my_data <- read.csv("data(1).csv")
my_data <- my_data[, c("US", "UK", "France", "Germany", "Italy")]

# --- DUMMY DATA FOR TESTING SCRIPT ---
set.seed(999)
dummy_data <- data.frame(
  US = cumsum(rnorm(500)),
  UK = cumsum(rnorm(500)),
  France = cumsum(rnorm(500)),
  Germany = cumsum(rnorm(500)),
  Italy = cumsum(rnorm(500))
)
# -------------------------------------

# Run Analysis
empirical_results <- run_country_analysis(dummy_data)

# Print cleanly to console
cat("\n=======================================================\n")
cat(" EMPIRICAL RESULTS (Replicating Tables 5 & 6) \n")
cat("=======================================================\n")
print(empirical_results)

# Optional: Save to CSV
# write.csv(empirical_results, "Empirical_DoubleCheck.csv", row.names = FALSE)