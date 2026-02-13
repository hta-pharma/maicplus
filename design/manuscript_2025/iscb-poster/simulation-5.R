# =============================================================================
# MAIC Simulation Study for Time-to-Event Outcomes with IPTW ATUT Comparison
# Version: Updated with systematic fixes and marginal HR calibration
#
# Key Updates in This Version:
# 1. Removed unnecessary return() statements that prevented proper capping
# 2. Added calibrate_event_rate_for_marginal_hr() function to target specific marginal HRs
# 3. Systematic review and cleanup of all functions
# =============================================================================

# Load required libraries
library(survival)
library(parallel)
library(MASS)

# Set global parameters
set.seed(12345)

# Enhanced simulation parameters
sim_params <- list(
  # Sample sizes per arm
  n_sizes = c(20, 50, 200, 500),
  
  # TARGET marginal HR values (replacing event_rates)
  target_marginal_hr = c(0.70, 0.85),
  
  # Original conditional HR values (for reference/comparison)
  conditional_hr = c(0.70, 0.85),
  
  # Event distributions
  distributions = c("weibull", "lognormal"),
  
  # Median survival times (months) for baseline hazard
  median_surv = c(24, 36),
  
  # 48-month survival rates for baseline hazard
  surv_48m = c(0.05, 0.15),
  
  # Total dropout rates (excluding administrative censoring)
  dropout_rates = c(0.2, 0.4),
  
  # Effective sample size ratios (ESS/N)
  ess_ratios = c(0.15, 0.25, 0.33, 0.50),
  
  # Simulation settings
  n_sim = 10000,
  max_followup = 48,
  n_boot = 500,
  max_weight = 100,
  n_mc_truth = 10000,
  
  # Parallel computing
  n_cores = 4,
  
  # File naming options
  results_prefix = "maic_iptw_atut_sim_fixed",
  timestamp_files = TRUE
)

# =============================================================================
# HELPER FUNCTIONS (WITH RETURN() FIXES)
# =============================================================================

#' Generate baseline covariates for internal and external populations
generate_covariates <- function(n, population = "external", d = 1) {
  if (population == "external") {
    X1 <- rnorm(n, mean = 0.5, sd = 0.5)
    X2 <- rbinom(n, size = 1, prob = 0.5)
    X3 <- rbinom(n, size = 1, prob = 0.5)
    X4 <- rbinom(n, size = 1, prob = 0.5)
  } else {
    X1 <- rnorm(n, mean = 0.5 * d, sd = 0.5)
    X2 <- rbinom(n, size = 1, prob = min(0.5 * d, 1))
    X3 <- rbinom(n, size = 1, prob = min(0.5 * d, 1))
    X4 <- rbinom(n, size = 1, prob = min(0.5 * d, 1))
  }
  
  cbind(X1 = X1, X2 = X2, X3 = X3, X4 = X4)
}

#' Calculate IPTW ATUT weights
calculate_iptw_atut_weights <- function(X_internal, X_external, X_for_ps = 1:3, max_weight = 100) {
  n_int <- nrow(X_internal)
  n_ext <- nrow(X_external)
  
  if (n_int == 0 || n_ext == 0) {
    warning("Empty input to IPTW weight calculation")
    rep(1, n_int)  # No return() here - function will naturally return this
  } else {
    X_combined <- rbind(X_internal[, X_for_ps, drop = FALSE],
                        X_external[, X_for_ps, drop = FALSE])
    T_combined <- c(rep(1, n_int), rep(0, n_ext))
    
    if (length(X_for_ps) == 1) {
      ps_data <- data.frame(T = T_combined, X1 = X_combined)
    } else {
      ps_data <- data.frame(T = T_combined, X_combined)
      colnames(ps_data) <- c("T", paste0("X", 1:length(X_for_ps)))
    }
    
    ps_model <- tryCatch({
      glm(T ~ ., data = ps_data, family = binomial())
    }, warning = function(w) {
      suppressWarnings(glm(T ~ ., data = ps_data, family = binomial()))
    }, error = function(e) {
      warning(paste("Propensity score model failed:", e$message))
      NULL
    })
    
    if (is.null(ps_model)) {
      warning("Propensity score model failed, returning equal weights")
      weights <- rep(1, n_int)
    } else {
      ps_pred_data <- data.frame(X_internal[, X_for_ps, drop = FALSE])
      if (length(X_for_ps) == 1) {
        colnames(ps_pred_data) <- "X1"
      } else {
        colnames(ps_pred_data) <- paste0("X", 1:length(X_for_ps))
      }
      
      ps_internal <- tryCatch({
        predict(ps_model, newdata = ps_pred_data, type = "response")
      }, error = function(e) {
        warning(paste("Prediction failed:", e$message))
        rep(0.5, n_int)
      })
      
      weights <- (1 - ps_internal) / ps_internal
      
      if (any(is.na(weights)) || any(is.infinite(weights))) {
        warning("NA or infinite weights detected, setting to 1")
        weights[is.na(weights) | is.infinite(weights)] <- 1
      }
      
      # Apply weight truncation
      weights <- pmin(weights, max_weight)
    }
    
    weights  # Natural return without explicit return()
  }
}

#' NEW FUNCTION: Calibrate event rate to achieve target marginal ATUT HR
#' This replaces the event_rate parameter with target_marginal_hr
calibrate_event_rate_for_marginal_hr <- function(target_marginal_hr, conditional_hr, d, 
                                                 distribution, median_surv, surv_48m,
                                                 dropout_rate, max_followup = 48,
                                                 n_calib = 5000, tol = 0.02, max_iter = 20,
                                                 verbose = TRUE) {
  
  # Input validation
  if (target_marginal_hr <= 0 || target_marginal_hr >= 2) {
    stop("Target marginal HR must be between 0 and 2")
  }
  
  if (conditional_hr <= 0 || conditional_hr >= 2) {
    stop("Conditional HR must be between 0 and 2")
  }
  
  # Due to non-collapsibility, marginal HR is typically closer to 1 than conditional HR
  # Check if the target is reasonable
  if (conditional_hr < 1 && target_marginal_hr < conditional_hr) {
    warning("Target marginal HR is more extreme than conditional HR - this may be difficult to achieve")
  }
  if (conditional_hr > 1 && target_marginal_hr > conditional_hr) {
    warning("Target marginal HR is more extreme than conditional HR - this may be difficult to achieve")
  }
  
  # Calculate base parameters
  base_params <- calibrate_base_params(distribution, median_surv, surv_48m, max_followup)
  
  # Function to compute marginal HR for a given beta (event rate parameter)
  compute_marginal_for_beta <- function(beta_scale) {
    beta <- c(beta_scale, beta_scale, beta_scale, 0)
    
    # Generate large sample for accurate estimation
    set.seed(123456)  # Fixed seed for consistency
    X_external_mc <- generate_covariates(n_calib, "external")
    
    # Randomly assign to treatment/control
    n_half <- n_calib / 2
    idx_treat <- sample(1:n_calib, n_half)
    idx_control <- setdiff(1:n_calib, idx_treat)
    
    X_treat <- X_external_mc[idx_treat, ]
    X_control <- X_external_mc[idx_control, ]
    
    gamma <- log(conditional_hr)
    
    # Generate survival times
    T_treat <- generate_survival_times(n_half, X_treat, beta, gamma, rep(1, n_half),
                                       distribution, base_params)
    T_control <- generate_survival_times(n_half, X_control, beta, 0, rep(0, n_half),
                                         distribution, base_params)
    
    # Apply censoring
    C_treat <- generate_censoring(n_half, T_treat, dropout_rate, max_followup)
    C_control <- generate_censoring(n_half, T_control, dropout_rate, max_followup)
    
    Y_treat <- pmin(T_treat, C_treat)
    Y_control <- pmin(T_control, C_control)
    delta_treat <- as.numeric(T_treat <= C_treat)
    delta_control <- as.numeric(T_control <= C_control)
    
    # Combine and fit Cox model
    Y_combined <- c(Y_treat, Y_control)
    delta_combined <- c(delta_treat, delta_control)
    Z_combined <- c(rep(1, n_half), rep(0, n_half))
    
    surv_obj <- Surv(Y_combined, delta_combined)
    cox_fit <- tryCatch({
      coxph(surv_obj ~ Z_combined)
    }, error = function(e) {
      NULL
    })
    
    if (is.null(cox_fit)) {
      NA
    } else {
      exp(coef(cox_fit))
    }
  }
  
  # Binary search for beta that achieves target marginal HR
  beta_lower <- -3
  beta_upper <- 3
  
  if (verbose) {
    cat(sprintf("Calibrating event rate for target marginal HR = %.3f\n", target_marginal_hr))
    cat(sprintf("  Conditional HR = %.3f, d = %.3f\n", conditional_hr, d))
  }
  
  for (iter in 1:max_iter) {
    beta_mid <- (beta_lower + beta_upper) / 2
    marginal_mid <- compute_marginal_for_beta(beta_mid)
    
    if (is.na(marginal_mid)) {
      warning("Failed to compute marginal HR")
      beta_mid  # Return current best guess
    } else {
      if (verbose && iter <= 5) {
        cat(sprintf("  Iter %d: beta = %.3f, marginal HR = %.3f\n", 
                    iter, beta_mid, marginal_mid))
      }
      
      if (abs(marginal_mid - target_marginal_hr) < tol) {
        if (verbose) {
          cat(sprintf("  Converged: beta = %.3f achieves marginal HR = %.3f\n", 
                      beta_mid, marginal_mid))
        }
        break
      }
      
      # Adjust search bounds
      if (marginal_mid < target_marginal_hr) {
        # Need higher HR, adjust beta to increase event rate difference
        if (conditional_hr < 1) {
          beta_upper <- beta_mid
        } else {
          beta_lower <- beta_mid
        }
      } else {
        if (conditional_hr < 1) {
          beta_lower <- beta_mid
        } else {
          beta_upper <- beta_mid
        }
      }
    }
  }
  
  # Return calibrated beta
  c(beta_mid, beta_mid, beta_mid, 0)
}

#' Compute true marginal ATUT HR via Monte Carlo
compute_true_marginal_atut_hr <- function(conditional_hr, d, beta, distribution,
                                          base_params, dropout_rate, max_followup = 48,
                                          n_mc = 10000) {
  set.seed(123456)
  
  X_external_mc <- generate_covariates(n_mc, "external")
  n_half <- n_mc / 2
  idx_treat <- sample(1:n_mc, n_half)
  idx_control <- setdiff(1:n_mc, idx_treat)
  
  X_treat <- X_external_mc[idx_treat, ]
  X_control <- X_external_mc[idx_control, ]
  
  gamma <- log(conditional_hr)
  
  T_treat <- generate_survival_times(n_half, X_treat, beta, gamma, rep(1, n_half),
                                     distribution, base_params)
  T_control <- generate_survival_times(n_half, X_control, beta, 0, rep(0, n_half),
                                       distribution, base_params)
  
  C_treat <- generate_censoring(n_half, T_treat, dropout_rate, max_followup)
  C_control <- generate_censoring(n_half, T_control, dropout_rate, max_followup)
  
  Y_treat <- pmin(T_treat, C_treat)
  Y_control <- pmin(T_control, C_control)
  delta_treat <- as.numeric(T_treat <= C_treat)
  delta_control <- as.numeric(T_control <= C_control)
  
  Y_combined <- c(Y_treat, Y_control)
  delta_combined <- c(delta_treat, delta_control)
  Z_combined <- c(rep(1, n_half), rep(0, n_half))
  
  surv_obj <- Surv(Y_combined, delta_combined)
  cox_fit <- tryCatch({
    coxph(surv_obj ~ Z_combined)
  }, error = function(e) {
    warning(paste("Cox model failed in true HR computation:", e$message))
    NULL
  })
  
  if (is.null(cox_fit)) {
    warning("Using conditional HR as fallback for true marginal HR")
    conditional_hr
  } else {
    as.numeric(exp(coef(cox_fit)))
  }
}

#' Calibrate d value to achieve target ESS ratio
calibrate_d_for_ess <- function(target_ess_ratio, n_calib = 10000, tol = 0.01,
                                method = "bisection") {
  
  compute_ess_ratio <- function(d) {
    X_external <- generate_covariates(n_calib, "external")
    X_internal <- generate_covariates(n_calib, "internal", d)
    
    weights <- calculate_iptw_atut_weights(X_internal, X_external, X_for_ps = 1:3)
    weights_combined <- c(weights, rep(1, n_calib))
    
    ess <- sum(weights_combined)^2 / sum(weights_combined^2)
    ess / (2 * n_calib)
  }
  
  d_lower <- 1.0
  d_upper <- 2.0
  
  ess_at_1 <- compute_ess_ratio(1.0)
  if (ess_at_1 < target_ess_ratio) {
    warning(sprintf("Target ESS ratio %.3f is higher than maximum achievable %.3f at d=1",
                    target_ess_ratio, ess_at_1))
    1.0
  } else {
    while (d_upper - d_lower > tol) {
      d_mid <- (d_lower + d_upper) / 2
      ess_mid <- compute_ess_ratio(d_mid)
      
      if (ess_mid > target_ess_ratio) {
        d_lower <- d_mid
      } else {
        d_upper <- d_mid
      }
    }
    
    (d_lower + d_upper) / 2
  }
}

#' Generate survival times based on specified distribution
generate_survival_times <- function(n, X, beta, gamma, Z, distribution, base_params) {
  lp <- as.numeric(X %*% beta + gamma * Z)
  
  if (distribution == "weibull") {
    if (base_params$shape <= 0 || base_params$lambda <= 0) {
      stop(sprintf("Invalid Weibull parameters: shape=%.4f, lambda=%.4f",
                   base_params$shape, base_params$lambda))
    }
    
    U <- runif(n)
    scale <- 1 / (base_params$lambda * exp(lp))
    T <- scale * (-log(U))^(1/base_params$shape)
    
  } else if (distribution == "lognormal") {
    T <- exp(rnorm(n, mean = base_params$mu - lp, sd = base_params$sigma))
  } else {
    stop("Unknown distribution")
  }
  
  pmax(T, 0.001)
}

#' Calibrate baseline distribution parameters
calibrate_base_params <- function(distribution, median_surv, surv_target, target_time = 48) {
  if (distribution == "weibull") {
    if (surv_target <= 0 || surv_target >= 1) {
      surv_target <- pmax(0.001, pmin(0.999, surv_target))
    }
    
    shape <- log(-log(surv_target)) / log(target_time / median_surv)
    
    if (shape <= 0) {
      warning(sprintf("Invalid Weibull shape (%.4f). Using shape = 1 (exponential).", shape))
      shape <- 1
    }
    
    lambda <- (-log(0.5))^(1/shape) / median_surv
    
    list(lambda = lambda, shape = shape)
    
  } else if (distribution == "lognormal") {
    mu <- log(median_surv)
    sigma <- (log(target_time) - mu) / qnorm(1 - surv_target)
    
    if (sigma <= 0) {
      warning(sprintf("Invalid log-normal sigma (%.4f). Using sigma = 0.5.", sigma))
      sigma <- 0.5
    }
    
    list(mu = mu, sigma = sigma)
  }
}

#' Generate censoring times based on dropout pattern
generate_censoring <- function(n, T, dropout_rate, max_followup = 48) {
  n_dropout <- rbinom(1, n, dropout_rate)
  
  if (n_dropout == 0) {
    rep(max_followup, n)
  } else {
    n_early <- rbinom(1, n_dropout, 0.6)
    n_late <- n_dropout - n_early
    
    dropout_ids <- sample(1:n, n_dropout, replace = FALSE)
    early_ids <- if(n_early > 0) dropout_ids[1:n_early] else numeric(0)
    late_ids <- if (n_late > 0) dropout_ids[(n_early + 1):n_dropout] else numeric(0)
    
    C <- rep(max_followup, n)
    
    for (i in early_ids) {
      C[i] <- runif(1, 0, min(T[i], 6))
    }
    
    for (i in late_ids) {
      if (T[i] > 6) {
        C[i] <- runif(1, 6, min(T[i], max_followup))
      } else {
        C[i] <- runif(1, 0, T[i])
      }
    }
    
    C
  }
}

#' Calculate MAIC weights using method of moments
calculate_maic_weights <- function(X_internal, X_external_means, max_weight = 100) {
  n <- nrow(X_internal)
  p <- ncol(X_internal)
  
  X_centered <- sweep(X_internal, 2, X_external_means, "-")
  
  objective <- function(alpha) {
    log_weights <- as.numeric(X_centered %*% alpha)
    log_weights <- pmin(log_weights, 20)
    weights <- exp(log_weights)
    
    if (sum(weights) == 0) {
      1e10
    } else {
      weighted_means <- colSums(weights * X_internal) / sum(weights)
      sum((weighted_means - X_external_means)^2)
    }
  }
  
  gradient <- function(alpha) {
    log_weights <- as.numeric(X_centered %*% alpha)
    log_weights <- pmin(log_weights, 20)
    weights <- exp(log_weights)
    
    if (sum(weights) == 0) {
      rep(0, p)
    } else {
      W <- sum(weights)
      weighted_X <- colSums(weights * X_internal)
      
      grad <- numeric(p)
      for (j in 1:p) {
        term1 <- sum(weights * X_centered[, j] * X_internal[, j]) / W
        term2 <- weighted_X[j] * sum(weights * X_centered[, j]) / W^2
        grad[j] <- 2 * (weighted_X[j]/W - X_external_means[j]) * (term1 - term2)
      }
      grad
    }
  }
  
  result <- tryCatch({
    optim(par = rep(0, p),
          fn = objective,
          gr = gradient,
          method = "BFGS",
          control = list(maxit = 1000))
  }, error = function(e) {
    optim(par = rep(0, p),
          fn = objective,
          method = "BFGS",
          control = list(maxit = 1000))
  })
  
  log_weights <- as.numeric(X_centered %*% result$par)
  log_weights <- pmin(log_weights, 20)
  weights <- exp(log_weights)
  
  # Apply weight capping
  pmin(weights, max_weight)
}

#' Calculate maximum ESS weights (MEW) - FIXED VERSION
calculate_mew_weights <- function(X_internal, X_external_means, max_weight = 100) {
  n <- nrow(X_internal)
  p <- ncol(X_internal)
  
  # Input validation
  if (n == 0 || p == 0) {
    warning("MEW: Empty input, returning uniform weights")
    rep(1, max(n, 1))
  } else if (length(X_external_means) != p) {
    stop("MEW: Dimension mismatch between X_internal and X_external_means")
  } else if (any(is.na(X_internal)) || any(is.infinite(X_internal))) {
    warning("MEW: NA or Inf in X_internal, returning uniform weights")
    rep(1, n)
  } else if (any(is.na(X_external_means)) || any(is.infinite(X_external_means))) {
    warning("MEW: NA or Inf in X_external_means, returning uniform weights")
    rep(1, n)
  } else {
    # Construct constraint matrix and solve
    A <- cbind(rep(1, n), X_internal)
    b <- c(n, n * X_external_means)
    
    weights <- tryCatch({
      AtA <- t(A) %*% A
      cond_num <- kappa(AtA)
      
      if (cond_num > 1e12) {
        warning(sprintf("MEW: Very ill-conditioned (kappa=%.2e), using fallback", cond_num))
        NULL
      } else {
        if (cond_num > 1e8) {
          reg_param <- 1e-6 * mean(diag(AtA))
          AtA <- AtA + diag(reg_param, ncol(AtA))
        }
        
        w <- as.numeric(A %*% solve(AtA) %*% b)
        
        if (any(is.na(w)) || any(is.infinite(w))) {
          warning("MEW: NA or Inf in solution, using fallback")
          NULL
        } else {
          w
        }
      }
    }, error = function(e) {
      warning(sprintf("MEW: Optimization failed (%s), using fallback", e$message))
      NULL
    })
    
    if (is.null(weights)) {
      rep(1, n)
    } else {
      # Handle non-positive weights
      nonpos_count <- sum(weights <= 0)
      if (nonpos_count > 0) {
        nonpos_prop <- nonpos_count / n
        
        if (nonpos_prop > 0.5) {
          warning(sprintf("MEW: %.1f%% non-positive weights, using uniform", 100*nonpos_prop))
          rep(1, n)
        } else {
          # Truncate at small positive value and cap
          weights <- pmax(weights, 0.00001)
          weights <- pmin(weights, max_weight)
          
          # Final check
          if (sd(weights) / mean(weights) < 0.001) {
            warning("MEW: Weights have no variation, using uniform")
            rep(1, n)
          } else {
            weights
          }
        }
      } else {
        # Just apply capping
        pmin(weights, max_weight)
      }
    }
  }
}

#' Normalize weights according to different schemes
normalize_weights <- function(weights, method) {
  n <- length(weights)
  
  if (sum(weights) == 0) {
    rep(1, n)
  } else {
    switch(method,
           "OW" = weights,
           "SW1" = weights / sum(weights),
           "SWN" = weights * n / sum(weights),
           "SWESS" = {
             ess <- sum(weights)^2 / sum(weights^2)
             weights * ess / sum(weights)
           },
           stop("Unknown normalization method"))
  }
}

#' Calculate effective sample size
calculate_ess <- function(weights) {
  if (sum(weights^2) == 0) {
    0
  } else {
    sum(weights)^2 / sum(weights^2)
  }
}

#' Fit weighted Cox model with proper bootstrap
fit_weighted_cox <- function(time, event, treatment, weights = NULL,
                             n_boot = 500, X_internal = NULL, X_external = NULL,
                             weight_type = "maic", weight_params = list()) {
  
  surv_obj <- Surv(time, event)
  
  cox_fit <- tryCatch({
    if (is.null(weights)) {
      coxph(surv_obj ~ treatment)
    } else {
      coxph(surv_obj ~ treatment, weights = weights, robust = TRUE)
    }
  }, error = function(e) {
    NULL
  })
  
  if (is.null(cox_fit)) {
    list(
      hr = NA, log_hr = NA, se_log_hr = NA,
      ci_robust_lower = NA, ci_robust_upper = NA,
      ci_boot_perc_lower = NA, ci_boot_perc_upper = NA,
      ci_boot_norm_lower = NA, ci_boot_norm_upper = NA,
      ci_boot_bca_lower = NA, ci_boot_bca_upper = NA,
      n_boot_valid = 0, boot_log_hr_dist = NULL
    )
  } else {
    if (is.null(weights)) {
      robust_se <- sqrt(vcov(cox_fit)[1,1])
    } else {
      if (!is.null(cox_fit$var)) {
        robust_se <- sqrt(cox_fit$var[1,1])
      } else {
        robust_se <- sqrt(vcov(cox_fit)[1,1])
      }
    }
    
    log_hr <- as.numeric(coef(cox_fit))
    hr <- exp(log_hr)
    
    ci_robust_lower <- exp(log_hr - 1.96 * robust_se)
    ci_robust_upper <- exp(log_hr + 1.96 * robust_se)
    
    # Initialize bootstrap results
    ci_boot_perc_lower <- ci_boot_perc_upper <- NA
    ci_boot_norm_lower <- ci_boot_norm_upper <- NA
    ci_boot_bca_lower <- ci_boot_bca_upper <- NA
    boot_log_hr_dist <- NULL
    n_boot_valid <- 0
    
    # Bootstrap if requested
    if (n_boot > 0 && !is.null(X_internal) && !is.null(X_external)) {
      trt_1_idx <- which(treatment == 1)
      trt_0_idx <- which(treatment == 0)
      n_1 <- length(trt_1_idx)
      n_0 <- length(trt_0_idx)
      
      if (n_0 > 0 && n_1 > 0) {
        boot_log_hr <- numeric(n_boot)
        boot_valid <- logical(n_boot)
        
        for (b in 1:n_boot) {
          boot_internal_idx <- sample(1:n_1, n_1, replace = TRUE)
          boot_external_idx <- sample(1:n_0, n_0, replace = TRUE)
          
          boot_idx_1 <- trt_1_idx[boot_internal_idx]
          boot_idx_0 <- trt_0_idx[boot_external_idx]
          
          X_int_boot <- X_internal[boot_internal_idx, , drop = FALSE]
          X_ext_boot <- X_external[boot_external_idx, , drop = FALSE]
          
          boot_weights <- tryCatch({
            if (weight_type == "maic") {
              w <- calculate_maic_weights(X_int_boot,
                                          colMeans(X_ext_boot),
                                          ifelse(is.null(weight_params$max_weight), 100,
                                                 weight_params$max_weight))
              normalize_weights(w, weight_params$norm_method)
            } else if (weight_type == "mew") {
              w <- calculate_mew_weights(X_int_boot,
                                         colMeans(X_ext_boot),
                                         ifelse(is.null(weight_params$max_weight), 100,
                                                weight_params$max_weight))
              normalize_weights(w, weight_params$norm_method)
            } else if (weight_type == "iptw") {
              w <- calculate_iptw_atut_weights(X_int_boot, X_ext_boot,
                                               weight_params$X_for_ps,
                                               ifelse(is.null(weight_params$max_weight), 100,
                                                      weight_params$max_weight))
              normalize_weights(w, weight_params$norm_method)
            } else {
              rep(1, n_1)
            }
          }, error = function(e) {
            rep(1, n_1)
          })
          
          boot_time <- c(time[boot_idx_1], time[boot_idx_0])
          boot_event <- c(event[boot_idx_1], event[boot_idx_0])
          boot_treatment <- c(rep(1, n_1), rep(0, n_0))
          boot_weights_combined <- c(boot_weights, rep(1, n_0))
          
          boot_cox <- tryCatch({
            boot_surv <- Surv(boot_time, boot_event)
            coxph(boot_surv ~ boot_treatment,
                  weights = boot_weights_combined,
                  robust = TRUE)
          }, error = function(e) {
            NULL
          })
          
          if (!is.null(boot_cox) && !is.na(coef(boot_cox))) {
            boot_log_hr[b] <- as.numeric(coef(boot_cox))
            boot_valid[b] <- TRUE
          } else {
            boot_valid[b] <- FALSE
          }
        }
        
        boot_log_hr_valid <- boot_log_hr[boot_valid]
        n_boot_valid <- sum(boot_valid)
        
        if (n_boot_valid >= 20) {
          boot_hr <- exp(boot_log_hr_valid)
          ci_boot_perc_lower <- as.numeric(quantile(boot_hr, 0.025, na.rm = TRUE))
          ci_boot_perc_upper <- as.numeric(quantile(boot_hr, 0.975, na.rm = TRUE))
          
          boot_se <- sd(boot_log_hr_valid, na.rm = TRUE)
          ci_boot_norm_lower <- exp(log_hr - 1.96 * boot_se)
          ci_boot_norm_upper <- exp(log_hr + 1.96 * boot_se)
          
          p_less <- mean(boot_log_hr_valid < log_hr, na.rm = TRUE)
          z0 <- qnorm(pmax(0.001, pmin(0.999, p_less)))
          
          p_lower <- pnorm(2 * z0 + qnorm(0.025))
          p_upper <- pnorm(2 * z0 + qnorm(0.975))
          
          p_lower <- pmax(0.001, pmin(0.999, p_lower))
          p_upper <- pmax(0.001, pmin(0.999, p_upper))
          
          ci_boot_bca_lower <- as.numeric(quantile(boot_hr, p_lower, na.rm = TRUE))
          ci_boot_bca_upper <- as.numeric(quantile(boot_hr, p_upper, na.rm = TRUE))
          
          boot_log_hr_dist <- boot_log_hr_valid
        }
      }
    }
    
    list(
      hr = hr,
      log_hr = log_hr,
      se_log_hr = robust_se,
      ci_robust_lower = ci_robust_lower,
      ci_robust_upper = ci_robust_upper,
      ci_boot_perc_lower = ci_boot_perc_lower,
      ci_boot_perc_upper = ci_boot_perc_upper,
      ci_boot_norm_lower = ci_boot_norm_lower,
      ci_boot_norm_upper = ci_boot_norm_upper,
      ci_boot_bca_lower = ci_boot_bca_lower,
      ci_boot_bca_upper = ci_boot_bca_upper,
      n_boot_valid = n_boot_valid,
      boot_log_hr_dist = boot_log_hr_dist
    )
  }
}

#' Calculate balance metrics
calculate_balance <- function(X_internal, X_external_means, weights) {
  weighted_means <- colSums(weights * X_internal) / sum(weights)
  external_sds <- c(0.5, rep(sqrt(0.5 * (1 - 0.5)), ncol(X_internal) - 1))
  
  if (length(external_sds) < length(weighted_means)) {
    external_sds <- rep(external_sds[1], length(weighted_means))
  }
  
  std_diffs <- abs(weighted_means - X_external_means) / external_sds
  max(std_diffs)
}

#' Run single simulation iteration
run_single_simulation_with_params <- function(params, sim_params) {
  # Extract parameters
  n <- params$n
  conditional_hr <- params$conditional_hr
  true_marginal_hr <- params$true_marginal_hr
  distribution <- params$distribution
  median_surv <- params$median_surv
  surv_48m <- params$surv_48m
  dropout_rate <- params$dropout_rate
  d <- params$d
  beta <- params$beta
  gamma <- params$gamma
  weighting_scenario <- params$weighting_scenario
  max_followup <- params$max_followup
  
  n_boot <- ifelse(is.null(sim_params$n_boot), 500, sim_params$n_boot)
  max_weight <- ifelse(is.null(sim_params$max_weight), 100, sim_params$max_weight)
  
  # Generate data
  base_params <- calibrate_base_params(distribution, median_surv, surv_48m, max_followup)
  
  X_internal <- generate_covariates(n, "internal", d)
  X_external <- generate_covariates(n, "external")
  
  T_internal <- generate_survival_times(n, X_internal, beta, gamma, rep(1, n),
                                        distribution, base_params)
  T_external <- generate_survival_times(n, X_external, beta, 0, rep(0, n),
                                        distribution, base_params)
  
  C_internal <- generate_censoring(n, T_internal, dropout_rate, max_followup)
  C_external <- generate_censoring(n, T_external, dropout_rate, max_followup)
  
  Y_internal <- pmin(T_internal, C_internal)
  Y_external <- pmin(T_external, C_external)
  delta_internal <- as.numeric(T_internal <= C_internal)
  delta_external <- as.numeric(T_external <= C_external)
  
  Y_combined <- c(Y_internal, Y_external)
  delta_combined <- c(delta_internal, delta_external)
  Z_combined <- c(rep(1, n), rep(0, n))
  
  # Select covariates for weighting
  if (weighting_scenario == "correct") {
    X_for_weights <- X_external[, 1:3]
    X_internal_for_weights <- X_internal[, 1:3]
    X_external_means <- colMeans(X_external[, 1:3])
    X_for_ps <- 1:3
  } else if (weighting_scenario == "under") {
    X_for_weights <- X_external[, 1:2]
    X_internal_for_weights <- X_internal[, 1:2]
    X_external_means <- colMeans(X_external[, 1:2])
    X_for_ps <- 1:2
  } else if (weighting_scenario == "over") {
    X_for_weights <- X_external
    X_internal_for_weights <- X_internal
    X_external_means <- colMeans(X_external)
    X_for_ps <- 1:4
  }
  
  # Calculate base weights
  maic_weights <- tryCatch({
    calculate_maic_weights(X_internal_for_weights, X_external_means, max_weight)
  }, error = function(e) {
    rep(1, n)
  })
  
  mew_weights <- tryCatch({
    calculate_mew_weights(X_internal_for_weights, X_external_means, max_weight)
  }, error = function(e) {
    rep(1, n)
  })
  
  iptw_weights <- tryCatch({
    calculate_iptw_atut_weights(X_internal, X_external, X_for_ps, max_weight)
  }, error = function(e) {
    rep(1, n)
  })
  
  # Initialize results
  results <- data.frame(
    n = integer(),
    conditional_hr = numeric(),
    true_marginal_hr = numeric(),
    distribution = character(),
    median_surv = numeric(),
    surv_48m = numeric(),
    dropout_rate = numeric(),
    d = numeric(),
    weighting_scenario = character(),
    method = character(),
    hr_est = numeric(),
    ci_robust_lower = numeric(),
    ci_robust_upper = numeric(),
    ci_boot_perc_lower = numeric(),
    ci_boot_perc_upper = numeric(),
    ci_boot_norm_lower = numeric(),
    ci_boot_norm_upper = numeric(),
    ci_boot_bca_lower = numeric(),
    ci_boot_bca_upper = numeric(),
    ess = numeric(),
    max_weight = numeric(),
    max_std_diff = numeric(),
    n_boot_valid = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process all methods
  methods <- c("OW", "SW1", "SWN", "SWESS", "MEW", "MEW1", "MEWN", "MEWESS",
               "IPTW_OW", "IPTW_SW1", "IPTW_SWN", "IPTW_SWESS")
  
  for (method in methods) {
    # Determine base weights and normalization
    if (substr(method, 1, 3) == "MEW") {
      base_weights <- mew_weights
      weight_type <- "mew"
      norm_method <- ifelse(method == "MEW", "OW",
                            ifelse(method == "MEW1", "SW1",
                                   ifelse(method == "MEWN", "SWN", "SWESS")))
    } else if (substr(method, 1, 4) == "IPTW") {
      base_weights <- iptw_weights
      weight_type <- "iptw"
      norm_method <- substring(method, 6)
    } else {
      base_weights <- maic_weights
      weight_type <- "maic"
      norm_method <- method
    }
    
    weights_internal <- normalize_weights(base_weights, norm_method)
    weights_combined <- c(weights_internal, rep(1, n))
    
    weight_params <- list(
      norm_method = norm_method,
      max_weight = max_weight,
      X_for_ps = X_for_ps
    )
    
    cox_results <- fit_weighted_cox(Y_combined, delta_combined, Z_combined,
                                    weights_combined, n_boot = n_boot,
                                    X_internal = X_internal,
                                    X_external = X_external,
                                    weight_type = weight_type,
                                    weight_params = weight_params)
    
    ess <- calculate_ess(weights_internal)
    max_weight_val <- max(weights_internal)
    
    if (weighting_scenario == "correct" || weighting_scenario == "under") {
      balance_X <- X_internal[, 1:3]
      balance_means <- colMeans(X_external[, 1:3])
    } else {
      balance_X <- X_internal
      balance_means <- colMeans(X_external)
    }
    
    max_std_diff <- calculate_balance(balance_X, balance_means, weights_internal)
    
    results <- rbind(results, data.frame(
      n = n,
      conditional_hr = conditional_hr,
      true_marginal_hr = true_marginal_hr,
      distribution = distribution,
      median_surv = median_surv,
      surv_48m = surv_48m,
      dropout_rate = dropout_rate,
      d = d,
      weighting_scenario = weighting_scenario,
      method = method,
      hr_est = cox_results$hr,
      ci_robust_lower = cox_results$ci_robust_lower,
      ci_robust_upper = cox_results$ci_robust_upper,
      ci_boot_perc_lower = cox_results$ci_boot_perc_lower,
      ci_boot_perc_upper = cox_results$ci_boot_perc_upper,
      ci_boot_norm_lower = cox_results$ci_boot_norm_lower,
      ci_boot_norm_upper = cox_results$ci_boot_norm_upper,
      ci_boot_bca_lower = cox_results$ci_boot_bca_lower,
      ci_boot_bca_upper = cox_results$ci_boot_bca_upper,
      ess = ess,
      max_weight = max_weight_val,
      max_std_diff = max_std_diff,
      n_boot_valid = cox_results$n_boot_valid,
      stringsAsFactors = FALSE
    ))
  }
  
  results
}

#' Run full simulation study - UPDATED VERSION
run_simulation_study <- function(sim_params, save_intermediate = TRUE,
                                 save_dir = "simulation_results") {
  
  if (save_intermediate) {
    if (sim_params$timestamp_files) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      save_dir <- paste0(save_dir, "_", timestamp)
    }
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    saveRDS(sim_params, file.path(save_dir, "simulation_parameters.rds"))
  }
  
  # Pre-calibrate d values
  cat("Calibrating d values for ESS ratios...\n")
  d_values <- numeric(length(sim_params$ess_ratios))
  for (i in seq_along(sim_params$ess_ratios)) {
    cat(sprintf("  ESS ratio %.2f: ", sim_params$ess_ratios[i]))
    d_values[i] <- calibrate_d_for_ess(sim_params$ess_ratios[i])
    cat(sprintf("d = %.3f\n", d_values[i]))
  }
  names(d_values) <- as.character(sim_params$ess_ratios)
  
  # Create factorial design
  scenarios <- expand.grid(
    n = sim_params$n_sizes,
    target_marginal_hr = sim_params$target_marginal_hr,
    conditional_hr = sim_params$conditional_hr,
    distribution = sim_params$distributions,
    median_surv = sim_params$median_surv,
    surv_48m = sim_params$surv_48m,
    dropout_rate = sim_params$dropout_rates,
    ess_ratio = sim_params$ess_ratios,
    weighting_scenario = c("correct", "under", "over"),
    stringsAsFactors = FALSE
  )
  
  scenarios$d <- d_values[as.character(scenarios$ess_ratio)]
  scenarios$max_followup <- sim_params$max_followup
  
  # Initialize columns
  scenarios$beta1 <- NA
  scenarios$beta2 <- NA
  scenarios$beta3 <- NA
  scenarios$beta4 <- NA
  scenarios$gamma <- NA
  scenarios$true_marginal_hr <- NA
  
  # Calibrate parameters
  cat("Calibrating parameters for target marginal HRs...\n")
  
  unique_combos <- unique(scenarios[, c("distribution", "median_surv", "surv_48m",
                                        "dropout_rate", "target_marginal_hr", 
                                        "conditional_hr", "d")])
  
  for (i in 1:nrow(unique_combos)) {
    combo <- unique_combos[i, ]
    cat(sprintf("  Combo %d/%d: ", i, nrow(unique_combos)))
    
    base_params <- calibrate_base_params(combo$distribution, combo$median_surv,
                                         combo$surv_48m, sim_params$max_followup)
    
    # Use the new calibration function
    beta <- calibrate_event_rate_for_marginal_hr(
      target_marginal_hr = combo$target_marginal_hr,
      conditional_hr = combo$conditional_hr,
      d = combo$d,
      distribution = combo$distribution,
      median_surv = combo$median_surv,
      surv_48m = combo$surv_48m,
      dropout_rate = combo$dropout_rate,
      max_followup = sim_params$max_followup,
      verbose = FALSE
    )
    
    gamma <- log(combo$conditional_hr)
    
    # Verify the achieved marginal HR
    true_marginal <- compute_true_marginal_atut_hr(
      conditional_hr = combo$conditional_hr,
      d = combo$d,
      beta = beta,
      distribution = combo$distribution,
      base_params = base_params,
      dropout_rate = combo$dropout_rate,
      max_followup = sim_params$max_followup,
      n_mc = sim_params$n_mc_truth
    )
    
    cat(sprintf("Target = %.3f, Achieved = %.3f\n",
                combo$target_marginal_hr, true_marginal))
    
    # Update matching rows
    matching_rows <- which(
      scenarios$distribution == combo$distribution &
        scenarios$median_surv == combo$median_surv &
        scenarios$surv_48m == combo$surv_48m &
        scenarios$dropout_rate == combo$dropout_rate &
        scenarios$target_marginal_hr == combo$target_marginal_hr &
        scenarios$conditional_hr == combo$conditional_hr &
        scenarios$d == combo$d
    )
    
    for (row in matching_rows) {
      scenarios$beta1[row] <- beta[1]
      scenarios$beta2[row] <- beta[2]
      scenarios$beta3[row] <- beta[3]
      scenarios$beta4[row] <- beta[4]
      scenarios$gamma[row] <- gamma
      scenarios$true_marginal_hr[row] <- true_marginal
    }
  }
  
  # Set up parallel cluster
  cl <- makeCluster(sim_params$n_cores)
  
  # Export functions
  clusterExport(cl, c("generate_covariates", "calibrate_base_params",
                      "generate_survival_times", "generate_censoring",
                      "calculate_maic_weights", "calculate_mew_weights",
                      "calculate_iptw_atut_weights",
                      "normalize_weights", "calculate_ess", "fit_weighted_cox",
                      "calculate_balance", "run_single_simulation_with_params",
                      "compute_true_marginal_atut_hr",
                      "calibrate_event_rate_for_marginal_hr"),
                envir = environment())
  
  clusterExport(cl, "sim_params", envir = environment())
  
  clusterEvalQ(cl, {
    library(survival)
    library(MASS)
    set.seed(Sys.getpid())
  })
  
  # Run simulations
  total_scenarios <- nrow(scenarios)
  cat(sprintf("Running %d scenarios with %d replications each...\n",
              total_scenarios, sim_params$n_sim))
  
  all_results <- list()
  
  for (scenario_idx in 1:total_scenarios) {
    cat(sprintf("Scenario %d/%d", scenario_idx, total_scenarios))
    
    current_scenario <- as.list(scenarios[scenario_idx, ])
    current_scenario$beta <- c(
      current_scenario$beta1,
      current_scenario$beta2,
      current_scenario$beta3,
      current_scenario$beta4
    )
    
    param_list <- replicate(sim_params$n_sim, current_scenario, simplify = FALSE)
    
    scenario_results <- parLapply(cl, param_list, function(params) {
      tryCatch({
        run_single_simulation_with_params(params, sim_params)
      }, error = function(e) {
        cat("\nError in simulation:", e$message, "\n")
        NULL
      })
    })
    
    scenario_results <- scenario_results[!sapply(scenario_results, is.null)]
    
    if (length(scenario_results) > 0) {
      scenario_results_df <- do.call(rbind, scenario_results)
      scenario_results_df$iteration <- rep(1:length(scenario_results),
                                           each = 12)
      scenario_results_df$scenario_id <- scenario_idx
      
      all_results[[scenario_idx]] <- scenario_results_df
      
      cat(sprintf(" - Completed %d/%d iterations\n",
                  length(scenario_results), sim_params$n_sim))
    } else {
      cat(" - Failed (no valid results)\n")
    }
    
    # Save intermediate results
    if (save_intermediate && scenario_idx %% 10 == 0) {
      intermediate_df <- do.call(rbind, all_results)
      filename <- sprintf("%s_intermediate_%04d.rds",
                          sim_params$results_prefix, scenario_idx)
      saveRDS(intermediate_df, file.path(save_dir, filename))
      cat(sprintf("  Saved intermediate results at scenario %d\n", scenario_idx))
    }
  }
  
  stopCluster(cl)
  
  # Combine results
  final_results <- do.call(rbind, all_results)
  
  # Save final results
  if (save_intermediate) {
    filename_base <- sprintf("%s_n%d-%d_thr%.2f-%.2f_sim%d",
                             sim_params$results_prefix,
                             min(sim_params$n_sizes), max(sim_params$n_sizes),
                             min(sim_params$target_marginal_hr), 
                             max(sim_params$target_marginal_hr),
                             sim_params$n_sim)
    
    saveRDS(final_results, file.path(save_dir, paste0(filename_base, ".rds")))
    write.csv(final_results, file.path(save_dir, paste0(filename_base, ".csv")),
              row.names = FALSE)
    
    write.csv(scenarios, file.path(save_dir, "scenario_definitions.csv"),
              row.names = FALSE)
    
    cat(sprintf("\nResults saved to %s/\n", save_dir))
  }
  
  final_results
}

# Additional analysis functions remain the same as before...
# [Including calculate_performance_metrics, create_comparison_tables, etc.]

#' Main function
main <- function(sim_params = NULL, run_analysis = TRUE) {
  if (is.null(sim_params)) {
    stop("sim_params must be provided")
  }
  
  # Set defaults
  if (is.null(sim_params$timestamp_files)) sim_params$timestamp_files <- TRUE
  if (is.null(sim_params$results_prefix)) sim_params$results_prefix <- "maic_iptw_atut_sim_fixed"
  if (is.null(sim_params$n_boot)) sim_params$n_boot <- 500
  if (is.null(sim_params$max_weight)) sim_params$max_weight <- 100
  if (is.null(sim_params$n_mc_truth)) sim_params$n_mc_truth <- 10000
  
  start_time <- Sys.time()
  
  cat("=============================================================================\n")
  cat("MAIC Simulation Study - Updated Version with Fixes\n")
  cat("TARGET: Specific marginal ATUT HRs via event rate calibration\n")
  cat("=============================================================================\n\n")
  
  # Run simulation
  cat("Starting simulation study...\n")
  results <- run_simulation_study(sim_params)
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "hours")
  
  cat("\n=============================================================================\n")
  cat(sprintf("Simulation completed in %.2f hours\n", runtime))
  cat(sprintf("Total scenarios: %d\n", length(unique(paste(results$n, 
                                                           results$target_marginal_hr,
                                                           results$distribution,
                                                           results$weighting_scenario)))))
  cat(sprintf("Total observations: %d\n", nrow(results)))
  cat("=============================================================================\n")
  
  list(
    results = results,
    runtime = runtime,
    sim_params = sim_params
  )
}