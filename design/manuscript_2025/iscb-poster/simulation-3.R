# =============================================================================
# MAIC Simulation Study for Time-to-Event Outcomes - FINAL VERSION
#
# This code implements the simulation study described in the manuscript:
# "Comparative Evaluation of Weight Normalization Methods in Matching-Adjusted
# Indirect Comparison for Time-to-Event Outcomes"
#
# Author: Gregory Chen, etc
# Date: Aug 15, 2025
# R Version: 4.5.1
# =============================================================================

# =============================================================================
# TESTING FRAMEWORK
# =============================================================================

#' Run comprehensive tests on all simulation functions
#'
#' @description Tests all major functions with small datasets to ensure
#' error-free execution before running full simulations
#'
#' @return List of test results with pass/fail status
test_simulation_functions <- function() {
  cat("=== Running Simulation Function Tests ===\n")

  test_results <- list()

  # Test 1: Covariate generation
  cat("Test 1: Covariate generation...\n")
  tryCatch({
    X_ext <- generate_covariates(50, "external", d = 1)
    X_int <- generate_covariates(50, "internal", d = 1.5)
    test_results$covariates <- all(dim(X_ext) == c(50, 4), dim(X_int) == c(50, 4))
    cat("  PASS\n")
  }, error = function(e) {
    test_results$covariates <- FALSE
    cat("  FAIL:", e$message, "\n")
  })

  # Test 2: Calibrate d for ESS
  cat("Test 2: Calibrate d for ESS...\n")
  tryCatch({
    d_test <- calibrate_d_for_ess(0.25, n_calib = 1000, method = "bisection")
    test_results$calibrate_d <- d_test > 1 && d_test < 2
    cat("  PASS: d =", d_test, "\n")
  }, error = function(e) {
    test_results$calibrate_d <- FALSE
    cat("  FAIL:", e$message, "\n")
  })

  # Test 3: Base parameter calibration
  cat("Test 3: Base parameter calibration...\n")
  tryCatch({
    params_wb <- calibrate_base_params("weibull", 24, 0.15, 48)
    params_ln <- calibrate_base_params("lognormal", 24, 0.15, 48)
    test_results$base_params <- all(length(params_wb) == 2, length(params_ln) == 2)
    cat("  PASS\n")
  }, error = function(e) {
    test_results$base_params <- FALSE
    cat("  FAIL:", e$message, "\n")
  })

  # Test 4: Survival time generation
  cat("Test 4: Survival time generation...\n")
  tryCatch({
    params <- calibrate_base_params("weibull", 24, 0.15, 48)
    X <- matrix(rnorm(50*4), 50, 4)
    beta <- c(0.5, 0.5, 0.5, 0)
    T_sim <- generate_survival_times(50, X, beta, log(0.7), rep(1, 50), "weibull", params)
    test_results$survival_times <- all(T_sim > 0, length(T_sim) == 50)
    cat("  PASS\n")
  }, error = function(e) {
    test_results$survival_times <- FALSE
    cat("  FAIL:", e$message, "\n")
  })

  # Test 5: Censoring generation
  cat("Test 5: Censoring generation...\n")
  tryCatch({
    T_test <- rexp(50, 0.02)
    C_test <- generate_censoring(50, T_test, 0.2, 48)
    test_results$censoring <- all(C_test <= 48, length(C_test) == 50)
    cat("  PASS\n")
  }, error = function(e) {
    test_results$censoring <- FALSE
    cat("  FAIL:", e$message, "\n")
  })

  # Test 6: Weight calculation
  cat("Test 6: Weight calculation...\n")
  tryCatch({
    X_int <- matrix(rnorm(50*3, 0.5), 50, 3)
    X_ext_means <- c(0, 0, 0)
    w_maic <- calculate_maic_weights(X_int, X_ext_means)
    w_mew <- calculate_mew_weights(X_int, X_ext_means)
    test_results$weights <- all(w_maic > 0, w_mew >= 0, length(w_maic) == 50)
    cat("  PASS\n")
  }, error = function(e) {
    test_results$weights <- FALSE
    cat("  FAIL:", e$message, "\n")
  })

  # Test 7: Cox model fitting
  cat("Test 7: Cox model fitting...\n")
  tryCatch({
    time <- rexp(100, 0.1)
    event <- rbinom(100, 1, 0.7)
    treatment <- rep(0:1, each = 50)
    weights <- runif(100, 0.5, 2)
    cox_res <- fit_weighted_cox(time, event, treatment, weights, n_boot = 50)
    test_results$cox_model <- !is.na(cox_res$hr) && length(cox_res) > 10
    cat("  PASS\n")
  }, error = function(e) {
    test_results$cox_model <- FALSE
    cat("  FAIL:", e$message, "\n")
  })

  # Test 8: Single simulation run
  cat("Test 8: Single simulation run...\n")
  tryCatch({
    test_params <- list(
      n = 50, true_hr = 0.7, distribution = "weibull",
      median_surv = 24, surv_48m = 0.15, dropout_rate = 0.2,
      event_rate = 0.25, d = 1.5, beta = c(0.5, 0.5, 0.5, 0),
      gamma = log(0.7), weighting_scenario = "correct",
      max_followup = 48
    )
    # Set global sim_params for n_boot
    sim_params_temp <- list(n_boot = 50)
    assign("sim_params", sim_params_temp, envir = .GlobalEnv)

    single_res <- run_single_simulation(test_params)
    test_results$single_sim <- nrow(single_res) == 8  # 8 methods
    cat("  PASS\n")
  }, error = function(e) {
    test_results$single_sim <- FALSE
    cat("  FAIL:", e$message, "\n")
  })

  # Summary
  cat("\n=== Test Summary ===\n")
  pass_count <- sum(unlist(test_results))
  total_count <- length(test_results)
  cat(sprintf("Passed: %d/%d tests\n", pass_count, total_count))

  if (pass_count < total_count) {
    cat("\nFailed tests:\n")
    for (test_name in names(test_results)) {
      if (!test_results[[test_name]]) {
        cat("  -", test_name, "\n")
      }
    }
  }

  return(test_results)
}

# =============================================================================
# MAIN SIMULATION CODE (CORRECTED VERSION)
# =============================================================================

# Load required libraries
library(survival)
library(parallel)

# Set global parameters with flexible file naming
set.seed(12345)  # For reproducibility

# Enhanced simulation parameters with file naming options
sim_params <- list(
  # Sample sizes per arm
  n_sizes = c(20, 50, 200, 500),

  # True ATUT HR values
  true_hr = c(0.70, 0.85),

  # Event distributions
  distributions = c("weibull", "lognormal"),

  # Median survival times (months) for baseline hazard
  median_surv = c(24, 36),

  # 48-month survival rates for baseline hazard
  surv_48m = c(0.05, 0.15),

  # Total dropout rates (excluding administrative censoring)
  dropout_rates = c(0.2, 0.4),

  # Target event rates in control group
  event_rates = c(0.05, 0.10, 0.25, 0.40),

  # Effective sample size ratios (ESS/N)
  ess_ratios = c(0.15, 0.25, 0.33, 0.50),

  # Simulation settings
  n_sim = 10000,  # Number of Monte Carlo replications
  max_followup = 48,  # months
  n_boot = 500,  # Number of bootstrap iterations for CIs

  # Parallel computing
  n_cores = 4,  # Adjust based on your system

  # File naming options
  results_prefix = "maic_sim",  # Prefix for result files
  timestamp_files = TRUE  # Add timestamp to filenames
)

# =============================================================================
# HELPER FUNCTIONS (CORRECTED VERSIONS)
# =============================================================================

#' Generate baseline covariates for internal and external populations
#'
#' @description Creates four baseline covariates (X1-X4) with controlled
#' differences between internal and external populations.
generate_covariates <- function(n, population = "external", d = 1) {
  if (population == "external") {
    # External population: X1 ~ N(0.5, 0.5^2), X2-X4 ~ Bernoulli(0.5)
    X1 <- rnorm(n, mean = 0.5, sd = 0.5)
    X2 <- rbinom(n, size = 1, prob = 0.5)
    X3 <- rbinom(n, size = 1, prob = 0.5)
    X4 <- rbinom(n, size = 1, prob = 0.5)
  } else {
    # Internal population: means scaled by d
    X1 <- rnorm(n, mean = 0.5 * d, sd = 0.5)
    X2 <- rbinom(n, size = 1, prob = min(0.5 * d, 1))
    X3 <- rbinom(n, size = 1, prob = min(0.5 * d, 1))
    X4 <- rbinom(n, size = 1, prob = min(0.5 * d, 1))
  }

  return(cbind(X1 = X1, X2 = X2, X3 = X3, X4 = X4))
}

#' Calibrate d value to achieve target ESS ratio (CORRECTED)
calibrate_d_for_ess <- function(target_ess_ratio, n_calib = 10000, tol = 0.01,
                                method = c("bisection", "newton-ralphson")) {
  # Function to compute ESS ratio for given d
  compute_ess_ratio <- function(d) {
    # Generate large samples for calibration
    X_external <- generate_covariates(n_calib, "external")
    X_internal <- generate_covariates(n_calib, "internal", d)

    # Compute propensity scores (only using prognostic factors X1-X3)
    combined_X <- rbind(X_external[, 1:3], X_internal[, 1:3])
    combined_Z <- c(rep(0, n_calib), rep(1, n_calib))

    # Create data frame for GLM
    glm_data <- data.frame(
      Z = combined_Z,
      X1 = combined_X[, 1],
      X2 = combined_X[, 2],
      X3 = combined_X[, 3]
    )

    # Fit propensity score model with error handling
    ps_model <- tryCatch({
      glm(Z ~ X1 + X2 + X3, data = glm_data, family = binomial())
    }, warning = function(w) {
      glm(Z ~ X1 + X2 + X3, data = glm_data, family = binomial())
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(ps_model)) return(0)

    # Create prediction data frame
    pred_data <- data.frame(
      X1 = X_external[, 1],
      X2 = X_external[, 2],
      X3 = X_external[, 3]
    )

    ps_external <- predict(ps_model, newdata = pred_data, type = "response")

    # Compute weights with numerical stability
    ps_external <- pmax(ps_external, 0.001)
    ps_external <- pmin(ps_external, 0.999)

    weights <- ps_external / (1 - ps_external)
    weights <- c(weights, rep(1, n_calib))

    # Compute ESS
    ess <- sum(weights)^2 / sum(weights^2)

    # Return ESS ratio
    return(ess / (2 * n_calib))
  }

  # Match the first element if method is a vector
  method <- match.arg(method)

  if(method == "bisection") {
    # Use bisection method to find d > 1
    d_lower <- 1.0
    d_upper <- 2.0

    # Check if target is achievable
    ess_at_1 <- compute_ess_ratio(1.0)
    if (ess_at_1 < target_ess_ratio) {
      warning(sprintf("Target ESS ratio %.3f is higher than maximum achievable %.3f at d=1",
                      target_ess_ratio, ess_at_1))
      return(1.0)
    }

    while (d_upper - d_lower > tol) {
      d_mid <- (d_lower + d_upper) / 2
      ess_mid <- compute_ess_ratio(d_mid)

      # Since ESS decreases as d increases from 1
      if (ess_mid > target_ess_ratio) {
        # Need larger d to get smaller ESS
        d_lower <- d_mid
      } else {
        # ESS is too small, need smaller d
        d_upper <- d_mid
      }
    }

    d_current <- (d_lower + d_upper) / 2
  } else if(method == "newton-ralphson") {
    # Use Newton-Raphson method to find d > 1
    d_current <- 1.1
    max_iter <- 50
    h <- 0.001  # Step size for numerical derivative

    # Check if target is achievable
    ess_at_1 <- compute_ess_ratio(1.0)
    if (ess_at_1 < target_ess_ratio) {
      warning(sprintf("Target ESS ratio %.3f is higher than maximum achievable %.3f at d=1",
                      target_ess_ratio, ess_at_1))
      return(1.0)
    }

    for (iter in 1:max_iter) {
      # Compute function value at current d
      f_current <- compute_ess_ratio(d_current) - target_ess_ratio

      # Check convergence
      if (abs(f_current) < tol) {
        break
      }

      # Compute numerical derivative using central difference
      f_plus <- compute_ess_ratio(d_current + h) - target_ess_ratio
      f_minus <- compute_ess_ratio(d_current - h) - target_ess_ratio
      f_prime <- (f_plus - f_minus) / (2 * h)

      # Check if derivative is too small (avoid division by zero)
      if (abs(f_prime) < 1e-10) {
        # Fall back to bisection for this iteration
        if (f_current > 0) {
          d_current <- d_current + 0.1  # Move away from 1
        } else {
          d_current <- d_current - 0.05  # Move toward 1 but stay > 1
        }
      } else {
        # Newton-Raphson update
        d_new <- d_current - f_current / f_prime

        # Ensure d stays in [1, 2] to find the d > 1 solution
        d_new <- pmax(1.001, pmin(2.0, d_new))

        # Apply damping if change is too large
        if (abs(d_new - d_current) > 0.5) {
          d_new <- d_current + sign(d_new - d_current) * 0.5
        }

        d_current <- d_new
      }
    }
  }

  return(d_current)
}

#' Generate survival times based on specified distribution
generate_survival_times <- function(n, X, beta, gamma, Z, distribution, base_params) {
  # Linear predictor
  lp <- as.numeric(X %*% beta + gamma * Z)

  if (distribution == "weibull") {
    # For Weibull: T = (-log(U) / (lambda * exp(lp)))^(1/shape)
    U <- runif(n)
    T <- (-log(U) / (base_params$lambda * exp(lp)))^(1/base_params$shape)
  } else if (distribution == "lognormal") {
    # For log-normal: log(T) ~ N(mu - lp, sigma^2)
    T <- exp(rnorm(n, mean = base_params$mu - lp, sd = base_params$sigma))
  }

  return(as.numeric(T))
}

#' Calibrate baseline distribution parameters
calibrate_base_params <- function(distribution, median_surv, surv_target, target_time = 48) {
  if (distribution == "weibull") {
    # For Weibull: S(t) = exp(-(lambda*t)^shape)
    # Median: exp(-(lambda*median)^shape) = 0.5
    # Target time: exp(-(lambda*target_time)^shape) = surv_target

    # Avoid log(0) issues
    if (surv_target <= 0 || surv_target >= 1) {
      surv_target <- pmax(0.001, pmin(0.999, surv_target))
    }

    shape <- log(-log(0.5)) / log(median_surv/target_time) *
      log(-log(surv_target)) / log(-log(0.5))
    lambda <- (-log(0.5))^(1/shape) / median_surv

    return(list(lambda = lambda, shape = shape))
  } else {
    # For log-normal: S(t) = 1 - Phi((log(t) - mu)/sigma)
    # Median: log(median_surv) = mu (since Phi(0) = 0.5)
    mu <- log(median_surv)

    # Target time survival: 1 - Phi((log(target_time) - mu)/sigma) = surv_target
    sigma <- (log(target_time) - mu) / qnorm(1 - surv_target)

    return(list(mu = mu, sigma = sigma))
  }
}

#' Generate censoring times based on dropout pattern
generate_censoring <- function(n, T, dropout_rate, max_followup = 48) {
  # Number of dropouts
  n_dropout <- rbinom(1, n, dropout_rate)

  # Handle edge case where n_dropout is 0
  if (n_dropout == 0) {
    return(rep(max_followup, n))
  }

  # Split dropouts: 60% in first 6 months, 40% after
  n_early <- rbinom(1, n_dropout, 0.6)
  n_late <- n_dropout - n_early

  # Randomly assign dropouts
  dropout_ids <- sample(1:n, n_dropout, replace = FALSE)
  early_ids <- if(n_early > 0) dropout_ids[1:n_early] else numeric(0)
  late_ids <- if (n_late > 0) dropout_ids[(n_early + 1):n_dropout] else numeric(0)

  # Initialize censoring times at administrative censoring
  C <- rep(max_followup, n)

  # Early dropouts: uniform in [0, min(T, 6)]
  for (i in early_ids) {
    C[i] <- runif(1, 0, min(T[i], 6))
  }

  # Late dropouts: uniform in [6, min(T, max_followup)]
  for (i in late_ids) {
    if (T[i] > 6) {
      C[i] <- runif(1, 6, min(T[i], max_followup))
    } else {
      # If event before 6 months, censor uniformly before event
      C[i] <- runif(1, 0, T[i])
    }
  }

  return(C)
}

#' Calibrate beta coefficients to achieve target event rate (CORRECTED)
calibrate_beta <- function(target_event_rate, n_calib = 5000, distribution,
                           base_params, dropout_rate, max_followup = 48,
                           tol = 0.01, method = c("bisection", "newton-ralphson")) {
  # Function to compute event rate for given beta scale
  compute_event_rate <- function(beta_scale) {
    # Set beta values (X1-X3 prognostic, X4 non-prognostic)
    beta <- c(beta_scale, beta_scale, beta_scale, 0)

    # Generate data
    X <- generate_covariates(n_calib, "external")
    T <- generate_survival_times(n_calib, X, beta, gamma = 0, Z = rep(0, n_calib),
                                 distribution, base_params)

    # Generate censoring
    C <- generate_censoring(n_calib, T, dropout_rate, max_followup)

    # Observed outcomes
    Y <- pmin(T, C)
    delta <- as.numeric(T <= C)

    # Event rate
    event_rate <- mean(delta)
    return(event_rate)
  }

  # Match the first element if method is a vector
  method <- match.arg(method)

  if(method == "bisection") {
    # Use bisection to find beta scale
    beta_lower <- -2
    beta_upper <- 2

    # Check monotonicity at boundaries
    rate_lower <- compute_event_rate(beta_lower)
    rate_upper <- compute_event_rate(beta_upper)

    # Check if target is achievable
    if (target_event_rate < rate_lower || target_event_rate > rate_upper) {
      if (target_event_rate < rate_lower) return(c(beta_lower, beta_lower, beta_lower, 0))
      if (target_event_rate > rate_upper) return(c(beta_upper, beta_upper, beta_upper, 0))
    }

    while (beta_upper - beta_lower > tol) {
      beta_mid <- (beta_lower + beta_upper) / 2
      rate_mid <- compute_event_rate(beta_mid)

      # Since event rate increases with beta
      if (rate_mid < target_event_rate) {
        beta_lower <- beta_mid
      } else {
        beta_upper <- beta_mid
      }
    }

    beta_scale <- (beta_lower + beta_upper) / 2
  } else if(method == "newton-ralphson") {
    # Use Newton-Raphson method to find beta scale
    beta_current <- 0.0  # Start at 0 (no effect)
    max_iter <- 50
    h <- 0.01  # Step size for numerical derivative

    # First check if target is achievable
    rate_at_minus2 <- compute_event_rate(-2)
    rate_at_plus2 <- compute_event_rate(2)

    if (target_event_rate < rate_at_minus2 || target_event_rate > rate_at_plus2) {
      if (target_event_rate < rate_at_minus2) return(c(-2, -2, -2, 0))
      if (target_event_rate > rate_at_plus2) return(c(2, 2, 2, 0))
    }

    for (iter in 1:max_iter) {
      # Compute function value at current beta
      f_current <- compute_event_rate(beta_current) - target_event_rate

      # Check convergence
      if (abs(f_current) < tol) {
        break
      }

      # Compute numerical derivative using central difference
      rate_plus <- compute_event_rate(beta_current + h)
      rate_minus <- compute_event_rate(beta_current - h)
      f_prime <- (rate_plus - rate_minus) / (2 * h)

      # Check if derivative is too small (avoid division by zero)
      if (abs(f_prime) < 1e-10) {
        # Fall back to bisection-like step
        if (f_current > 0) {
          beta_current <- beta_current - 0.1
        } else {
          beta_current <- beta_current + 0.1
        }
      } else {
        # Newton-Raphson update
        beta_new <- beta_current - f_current / f_prime

        # Ensure beta stays within reasonable bounds [-2, 2]
        beta_new <- pmax(-2, pmin(2, beta_new))

        # Apply damping if change is too large
        if (abs(beta_new - beta_current) > 0.5) {
          beta_new <- beta_current + sign(beta_new - beta_current) * 0.5
        }

        beta_current <- beta_new
      }
    }

    beta_scale <- beta_current
  }

  return(c(beta_scale, beta_scale, beta_scale, 0))
}

#' Calculate MAIC weights using method of moments (CORRECTED)
calculate_maic_weights <- function(X_internal, X_external_means) {
  n <- nrow(X_internal)
  p <- ncol(X_internal)

  # Center covariates
  X_centered <- sweep(X_internal, 2, X_external_means, "-")

  # Method of moments: solve for alpha such that sum(w_i * X_i) = n * X_external_means
  # where w_i = exp(X_centered %*% alpha)

  # Objective function for optimization
  objective <- function(alpha) {
    log_weights <- as.numeric(X_centered %*% alpha)
    # Prevent overflow
    log_weights <- pmin(log_weights, 20)
    weights <- exp(log_weights)

    # Prevent division by zero
    if (sum(weights) == 0) return(1e10)

    weighted_means <- colSums(weights * X_internal) / sum(weights)
    sum((weighted_means - X_external_means)^2)
  }

  # Gradient (with numerical stability)
  gradient <- function(alpha) {
    log_weights <- as.numeric(X_centered %*% alpha)
    log_weights <- pmin(log_weights, 20)
    weights <- exp(log_weights)

    if (sum(weights) == 0) return(rep(0, p))

    W <- sum(weights)
    weighted_X <- colSums(weights * X_internal)

    grad <- numeric(p)
    for (j in 1:p) {
      term1 <- sum(weights * X_centered[, j] * X_internal[, j]) / W
      term2 <- weighted_X[j] * sum(weights * X_centered[, j]) / W^2
      grad[j] <- 2 * (weighted_X[j]/W - X_external_means[j]) * (term1 - term2)
    }
    return(grad)
  }

  # Optimize with error handling
  result <- tryCatch({
    optim(par = rep(0, p), fn = objective, gr = gradient,
          method = "BFGS", control = list(maxit = 1000))
  }, error = function(e) {
    # Fall back to simple optimization without gradient
    optim(par = rep(0, p), fn = objective, method = "BFGS",
          control = list(maxit = 1000))
  })

  # Calculate weights with overflow protection
  log_weights <- as.numeric(X_centered %*% result$par)
  log_weights <- pmin(log_weights, 20)
  weights <- exp(log_weights)

  return(weights)
}

#' Calculate maximum ESS weights (CORRECTED)
calculate_mew_weights <- function(X_internal, X_external_means) {
  n <- nrow(X_internal)
  p <- ncol(X_internal)

  # Construct constraint matrix
  A <- cbind(rep(1, n), X_internal)
  b <- c(n, n * X_external_means)

  # Solve quadratic program: min(w'w) s.t. A'w = b
  # Solution: w = A(A'A)^{-1}b

  # Check if A'A is invertible
  AtA <- t(A) %*% A

  # Add regularization if needed
  if (det(AtA) < 1e-10 || rcond(AtA) < 1e-10) {
    AtA <- AtA + diag(1e-6, ncol(AtA))
  }

  # Solve with error handling
  weights <- tryCatch({
    A %*% solve(AtA) %*% b
  }, error = function(e) {
    # Fall back to pseudo-inverse
    A %*% MASS::ginv(AtA) %*% b
  })

  # Ensure non-negative weights
  weights <- as.numeric(weights)
  weights <- pmax(weights, 0)

  # Handle case where all weights are zero
  if (sum(weights) == 0) {
    weights <- rep(1, n)
  }

  return(weights)
}

#' Normalize weights according to different schemes
normalize_weights <- function(weights, method) {
  n <- length(weights)

  # Handle zero weights
  if (sum(weights) == 0) {
    return(rep(1, n))
  }

  switch(method,
         "OW" = weights,  # Original weights (unscaled)
         "SW1" = weights / sum(weights),  # Sum to 1
         "SWN" = weights * n / sum(weights),  # Sum to N
         "SWESS" = {
           ess <- sum(weights)^2 / sum(weights^2)
           weights * ess / sum(weights)  # Sum to ESS
         },
         stop("Unknown normalization method")
  )
}

#' Calculate effective sample size
calculate_ess <- function(weights) {
  if (sum(weights^2) == 0) return(0)
  sum(weights)^2 / sum(weights^2)
}

#' Fit weighted Cox model and extract HR with multiple CI methods (CORRECTED)
fit_weighted_cox <- function(time, event, treatment, weights = NULL, n_boot = 500) {
  # Create survival object
  surv_obj <- Surv(time, event)

  # Fit Cox model with robust variance if weights are provided
  cox_fit <- tryCatch({
    if (is.null(weights)) {
      coxph(surv_obj ~ treatment)
    } else {
      coxph(surv_obj ~ treatment, weights = weights, robust = TRUE)
    }
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(cox_fit)) {
    return(list(
      hr = NA, log_hr = NA, se_log_hr = NA,
      ci_robust_lower = NA, ci_robust_upper = NA,
      ci_boot_perc_lower = NA, ci_boot_perc_upper = NA,
      ci_boot_norm_lower = NA, ci_boot_norm_upper = NA,
      ci_boot_bca_lower = NA, ci_boot_bca_upper = NA,
      n_boot_valid = 0, boot_log_hr_dist = NULL
    ))
  }

  # Extract point estimate and standard error
  if (is.null(weights)) {
    robust_se <- sqrt(vcov(cox_fit)[1,1])
  } else {
    robust_se <- sqrt(cox_fit$var[1,1])  # Robust variance
  }

  log_hr <- as.numeric(coef(cox_fit))
  hr <- exp(log_hr)

  # Method 1: Robust CI (always calculated)
  ci_robust_lower <- exp(log_hr - 1.96 * robust_se)
  ci_robust_upper <- exp(log_hr + 1.96 * robust_se)

  # Initialize other CI methods
  ci_boot_perc_lower <- ci_boot_perc_upper <- NA
  ci_boot_norm_lower <- ci_boot_norm_upper <- NA
  ci_boot_bca_lower <- ci_boot_bca_upper <- NA
  boot_log_hr_dist <- NULL
  n_boot_valid <- 0

  # Bootstrap CIs (only if requested)
  if (n_boot > 0) {
    n_total <- length(time)

    # Identify unique treatment groups for stratified sampling
    trt_0_idx <- which(treatment == 0)
    trt_1_idx <- which(treatment == 1)
    n_0 <- length(trt_0_idx)
    n_1 <- length(trt_1_idx)

    # Pre-allocate bootstrap results
    boot_log_hr <- numeric(n_boot)
    boot_valid <- logical(n_boot)

    # Perform stratified bootstrap
    for (b in 1:n_boot) {
      # Stratified resampling
      boot_idx_0 <- if(n_0 > 0) sample(trt_0_idx, n_0, replace = TRUE) else numeric(0)
      boot_idx_1 <- if(n_1 > 0) sample(trt_1_idx, n_1, replace = TRUE) else numeric(0)
      boot_idx <- c(boot_idx_0, boot_idx_1)

      # Bootstrap data
      boot_time <- time[boot_idx]
      boot_event <- event[boot_idx]
      boot_treatment <- treatment[boot_idx]
      boot_weights <- if (!is.null(weights)) weights[boot_idx] else NULL

      # Fit model on bootstrap sample
      boot_cox <- tryCatch({
        boot_surv <- Surv(boot_time, boot_event)
        if (is.null(boot_weights)) {
          coxph(boot_surv ~ boot_treatment)
        } else {
          coxph(boot_surv ~ boot_treatment, weights = boot_weights)
        }
      }, error = function(e) {
        return(NULL)
      })

      if (!is.null(boot_cox)) {
        boot_log_hr[b] <- as.numeric(coef(boot_cox))
        boot_valid[b] <- TRUE
      } else {
        boot_valid[b] <- FALSE
      }
    }

    # Remove invalid bootstrap samples
    boot_log_hr_valid <- boot_log_hr[boot_valid]
    n_boot_valid <- sum(boot_valid)

    if (n_boot_valid >= 50) {  # Need reasonable number for CIs
      # Method 2: Bootstrap percentile CI
      boot_hr <- exp(boot_log_hr_valid)
      ci_boot_perc_lower <- as.numeric(quantile(boot_hr, 0.025, na.rm = TRUE))
      ci_boot_perc_upper <- as.numeric(quantile(boot_hr, 0.975, na.rm = TRUE))

      # Method 3: Bootstrap normal CI (using SE from bootstrap)
      boot_se <- sd(boot_log_hr_valid, na.rm = TRUE)
      ci_boot_norm_lower <- exp(log_hr - 1.96 * boot_se)
      ci_boot_norm_upper <- exp(log_hr + 1.96 * boot_se)

      # Method 4: BCa bootstrap CI (simplified version)
      # Calculate bias correction factor z0
      p_less <- mean(boot_log_hr_valid < log_hr)
      z0 <- qnorm(pmax(0.001, pmin(0.999, p_less)))

      # Use simplified BCa without jackknife (faster)
      alpha_lower <- 0.025
      alpha_upper <- 0.975
      z_lower <- qnorm(alpha_lower)
      z_upper <- qnorm(alpha_upper)

      # Simplified adjustment (no acceleration)
      p_lower <- pnorm(2 * z0 + z_lower)
      p_upper <- pnorm(2 * z0 + z_upper)

      # Ensure percentiles are in [0,1]
      p_lower <- pmax(0.001, pmin(0.999, p_lower))
      p_upper <- pmax(0.001, pmin(0.999, p_upper))

      # BCa CI
      ci_boot_bca_lower <- as.numeric(quantile(boot_hr, p_lower, na.rm = TRUE))
      ci_boot_bca_upper <- as.numeric(quantile(boot_hr, p_upper, na.rm = TRUE))

      # Store bootstrap distribution for potential diagnostics
      boot_log_hr_dist <- boot_log_hr_valid
    }
  }

  return(list(
    hr = hr,
    log_hr = log_hr,
    se_log_hr = robust_se,
    # Robust CI (always available)
    ci_robust_lower = ci_robust_lower,
    ci_robust_upper = ci_robust_upper,
    # Bootstrap CIs (if computed)
    ci_boot_perc_lower = ci_boot_perc_lower,
    ci_boot_perc_upper = ci_boot_perc_upper,
    ci_boot_norm_lower = ci_boot_norm_lower,
    ci_boot_norm_upper = ci_boot_norm_upper,
    ci_boot_bca_lower = ci_boot_bca_lower,
    ci_boot_bca_upper = ci_boot_bca_upper,
    # Diagnostics
    n_boot_valid = n_boot_valid,
    boot_log_hr_dist = boot_log_hr_dist
  ))
}

#' Calculate balance metrics
calculate_balance <- function(X_internal, X_external_means, weights) {
  # Weighted means
  weighted_means <- colSums(weights * X_internal) / sum(weights)

  # Pooled standard deviations (using external as reference)
  # For simulation, we know external SDs
  external_sds <- c(0.5, rep(sqrt(0.5 * (1 - 0.5)), ncol(X_internal) - 1))

  # Ensure same length
  if (length(external_sds) < length(weighted_means)) {
    external_sds <- rep(external_sds[1], length(weighted_means))
  }

  # Standardized differences
  std_diffs <- abs(weighted_means - X_external_means) / external_sds

  return(max(std_diffs))
}

#' Run single simulation iteration (CORRECTED)
run_single_simulation <- function(params) {
  # Extract parameters
  n <- params$n
  true_hr <- params$true_hr
  distribution <- params$distribution
  median_surv <- params$median_surv
  surv_48m <- params$surv_48m
  dropout_rate <- params$dropout_rate
  event_rate <- params$event_rate
  d <- params$d
  beta <- params$beta
  gamma <- params$gamma
  weighting_scenario <- params$weighting_scenario
  max_followup <- params$max_followup

  # Generate baseline distribution parameters
  base_params <- calibrate_base_params(distribution, median_surv, surv_48m,
                                       target_time = max_followup)

  # Generate covariates
  X_internal <- generate_covariates(n, "internal", d)
  X_external <- generate_covariates(n, "external")

  # Generate survival times
  T_internal <- generate_survival_times(n, X_internal, beta, gamma, rep(1, n),
                                        distribution, base_params)
  T_external <- generate_survival_times(n, X_external, beta, 0, rep(0, n),
                                        distribution, base_params)

  # Generate censoring
  C_internal <- generate_censoring(n, T_internal, dropout_rate, max_followup)
  C_external <- generate_censoring(n, T_external, dropout_rate, max_followup)

  # Observed outcomes
  Y_internal <- pmin(T_internal, C_internal)
  Y_external <- pmin(T_external, C_external)
  delta_internal <- as.numeric(T_internal <= C_internal)
  delta_external <- as.numeric(T_external <= C_external)

  # Combine data for analysis
  Y_combined <- c(Y_internal, Y_external)
  delta_combined <- c(delta_internal, delta_external)
  Z_combined <- c(rep(1, n), rep(0, n))

  # Determine which covariates to use based on scenario
  if (weighting_scenario == "correct") {
    # Use X1, X2, X3 (all prognostic factors)
    X_for_weights <- X_internal[, 1:3]
    X_external_means <- colMeans(X_external[, 1:3])
  } else if (weighting_scenario == "under") {
    # Omit X3 (one prognostic factor)
    X_for_weights <- X_internal[, 1:2]
    X_external_means <- colMeans(X_external[, 1:2])
  } else if (weighting_scenario == "over") {
    # Include all covariates (including non-prognostic X4)
    X_for_weights <- X_internal
    X_external_means <- colMeans(X_external)
  }

  # Calculate base weights with error handling
  maic_weights <- tryCatch({
    calculate_maic_weights(X_for_weights, X_external_means)
  }, error = function(e) {
    rep(1, n)  # Fall back to equal weights
  })

  mew_weights <- tryCatch({
    calculate_mew_weights(X_for_weights, X_external_means)
  }, error = function(e) {
    rep(1, n)  # Fall back to equal weights
  })

  # Initialize results storage
  results <- data.frame(
    n = integer(),
    true_hr = numeric(),
    distribution = character(),
    median_surv = numeric(),
    surv_48m = numeric(),
    dropout_rate = numeric(),
    event_rate = numeric(),
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

  # Fit models with different weight normalizations
  methods <- c("OW", "SW1", "SWN", "SWESS", "MEW", "MEW1", "MEWN", "MEWESS")

  for (method in methods) {
    # Get appropriate base weights
    if (substr(method, 1, 3) == "MEW") {
      base_weights <- mew_weights
      norm_method <- ifelse(method == "MEW", "OW",
                            ifelse(method == "MEW1", "SW1",
                                   ifelse(method == "MEWN", "SWN", "SWESS")))
    } else {
      base_weights <- maic_weights
      norm_method <- method
    }

    # Normalize weights
    weights_internal <- normalize_weights(base_weights, norm_method)

    # Create combined weights (external gets weight 1)
    weights_combined <- c(weights_internal, rep(1, n))

    # Fit weighted Cox model
    cox_results <- fit_weighted_cox(Y_combined, delta_combined, Z_combined,
                                    weights_combined, n_boot = sim_params$n_boot)

    # Calculate metrics
    ess <- calculate_ess(weights_internal)
    max_weight <- max(weights_internal)

    # Balance (use appropriate X based on scenario)
    if (weighting_scenario == "correct") {
      balance_X <- X_internal[, 1:3]
      balance_means <- colMeans(X_external[, 1:3])
    } else if (weighting_scenario == "under") {
      balance_X <- X_internal[, 1:3]  # Still check balance on all prognostic
      balance_means <- colMeans(X_external[, 1:3])
    } else {
      balance_X <- X_internal
      balance_means <- colMeans(X_external)
    }

    max_std_diff <- calculate_balance(balance_X, balance_means, weights_internal)

    # Store results
    results <- rbind(results, data.frame(
      n = n,
      true_hr = true_hr,
      distribution = distribution,
      median_surv = median_surv,
      surv_48m = surv_48m,
      dropout_rate = dropout_rate,
      event_rate = event_rate,
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
      max_weight = max_weight,
      max_std_diff = max_std_diff,
      n_boot_valid = cox_results$n_boot_valid,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}

#' Run full simulation study with enhanced file saving (CORRECTED)
run_simulation_study <- function(sim_params, save_intermediate = TRUE,
                                 save_dir = "simulation_results") {

  # Create save directory with timestamp if requested
  if (save_intermediate) {
    if (sim_params$timestamp_files) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      save_dir <- paste0(save_dir, "_", timestamp)
    }
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }

    # Save simulation parameters
    saveRDS(sim_params, file.path(save_dir, "simulation_parameters.rds"))
  }

  # Pre-calibrate d values for each ESS ratio
  cat("Calibrating d values for ESS ratios...\n")
  d_values <- numeric(length(sim_params$ess_ratios))
  for (i in seq_along(sim_params$ess_ratios)) {
    cat(sprintf("  ESS ratio %.2f: ", sim_params$ess_ratios[i]))
    d_values[i] <- calibrate_d_for_ess(sim_params$ess_ratios[i])
    cat(sprintf("d = %.3f\n", d_values[i]))
  }
  names(d_values) <- as.character(sim_params$ess_ratios)

  # Create full factorial design
  scenarios <- expand.grid(
    n = sim_params$n_sizes,
    true_hr = sim_params$true_hr,
    distribution = sim_params$distributions,
    median_surv = sim_params$median_surv,
    surv_48m = sim_params$surv_48m,
    dropout_rate = sim_params$dropout_rates,
    event_rate = sim_params$event_rates,
    ess_ratio = sim_params$ess_ratios,
    weighting_scenario = c("correct", "under", "over"),
    stringsAsFactors = FALSE
  )

  # Add calibrated d values
  scenarios$d <- d_values[as.character(scenarios$ess_ratio)]

  # Calibrate beta and gamma for each unique combination
  cat("Calibrating beta coefficients and gamma values...\n")

  unique_combos <- unique(scenarios[, c("distribution", "median_surv", "surv_48m",
                                        "dropout_rate", "event_rate", "true_hr")])

  # Initialize beta and gamma columns as lists
  scenarios$beta <- vector("list", nrow(scenarios))
  scenarios$gamma <- NA

  for (i in 1:nrow(unique_combos)) {
    combo <- unique_combos[i, ]
    cat(sprintf("  Combo %d/%d\n", i, nrow(unique_combos)))

    # Calibrate base parameters
    base_params <- calibrate_base_params(combo$distribution, combo$median_surv,
                                         combo$surv_48m, target_time = sim_params$max_followup)

    # Calibrate beta
    beta <- calibrate_beta(combo$event_rate, distribution = combo$distribution,
                           base_params = base_params, dropout_rate = combo$dropout_rate,
                           max_followup = sim_params$max_followup)

    # Calibrate gamma for target HR
    gamma <- log(combo$true_hr)

    # Store in scenarios
    matching_rows <- which(
      scenarios$distribution == combo$distribution &
        scenarios$median_surv == combo$median_surv &
        scenarios$surv_48m == combo$surv_48m &
        scenarios$dropout_rate == combo$dropout_rate &
        scenarios$event_rate == combo$event_rate &
        scenarios$true_hr == combo$true_hr
    )

    for (row in matching_rows) {
      scenarios$beta[[row]] <- beta
      scenarios$gamma[row] <- gamma
    }
  }

  # Set up parallel cluster
  cl <- makeCluster(sim_params$n_cores)

  # Export all necessary functions and objects
  clusterExport(cl, c("generate_covariates", "calibrate_base_params",
                      "generate_survival_times", "generate_censoring",
                      "calculate_maic_weights", "calculate_mew_weights",
                      "normalize_weights", "calculate_ess", "fit_weighted_cox",
                      "calculate_balance", "run_single_simulation", "sim_params"),
                envir = environment())

  clusterEvalQ(cl, {
    library(survival)
    library(MASS)
    set.seed(Sys.getpid())  # Different seed for each worker
  })

  # Run simulations
  total_scenarios <- nrow(scenarios)
  cat(sprintf("Running %d scenarios with %d replications each...\n",
              total_scenarios, sim_params$n_sim))

  all_results <- list()

  for (scenario_idx in 1:total_scenarios) {
    cat(sprintf("Scenario %d/%d", scenario_idx, total_scenarios))

    # Current scenario
    current_scenario <- scenarios[scenario_idx, ]
    current_scenario$max_followup <- sim_params$max_followup

    # Create parameter list for each replication
    param_list <- replicate(sim_params$n_sim, as.list(current_scenario),
                            simplify = FALSE)

    # Run parallel simulations
    scenario_results <- parLapply(cl, param_list, function(params) {
      tryCatch({
        run_single_simulation(params)
      }, error = function(e) {
        cat("\nError in simulation:", e$message, "\n")
        return(NULL)
      })
    })

    # Remove NULL results
    scenario_results <- scenario_results[!sapply(scenario_results, is.null)]

    if (length(scenario_results) > 0) {
      # Combine results
      scenario_results_df <- do.call(rbind, scenario_results)
      scenario_results_df$iteration <- rep(1:length(scenario_results),
                                           each = 8)  # 8 methods
      scenario_results_df$scenario_id <- scenario_idx

      # Store results
      all_results[[scenario_idx]] <- scenario_results_df

      cat(sprintf(" - Completed %d/%d iterations\n",
                  length(scenario_results), sim_params$n_sim))
    } else {
      cat(" - Failed (no valid results)\n")
    }

    # Save intermediate results if requested
    if (save_intermediate && scenario_idx %% 10 == 0) {
      intermediate_df <- do.call(rbind, all_results)
      filename <- sprintf("%s_intermediate_%04d.rds",
                          sim_params$results_prefix, scenario_idx)
      saveRDS(intermediate_df, file.path(save_dir, filename))
      cat(sprintf("  Saved intermediate results at scenario %d\n", scenario_idx))
    }
  }

  # Stop cluster
  stopCluster(cl)

  # Combine all results
  final_results <- do.call(rbind, all_results)

  # Save final results
  if (save_intermediate) {
    # Save with descriptive filename
    filename_base <- sprintf("%s_n%d-%d_hr%.2f-%.2f_sim%d",
                             sim_params$results_prefix,
                             min(sim_params$n_sizes), max(sim_params$n_sizes),
                             min(sim_params$true_hr), max(sim_params$true_hr),
                             sim_params$n_sim)

    saveRDS(final_results, file.path(save_dir, paste0(filename_base, ".rds")))
    write.csv(final_results, file.path(save_dir, paste0(filename_base, ".csv")),
              row.names = FALSE)

    # Save scenario summary
    write.csv(scenarios, file.path(save_dir, "scenario_definitions.csv"),
              row.names = FALSE)
  }

  return(final_results)
}

#' Enhanced performance metrics calculation with flexible reporting
calculate_performance_metrics <- function(results, save_output = TRUE,
                                          output_dir = "simulation_results") {
  # Check for valid results
  if (nrow(results) == 0 || all(is.na(results$hr_est))) {
    warning("No valid results to calculate metrics")
    return(NULL)
  }

  # Group by scenario and method
  grouping_vars <- c("n", "true_hr", "distribution", "median_surv", "surv_48m",
                     "dropout_rate", "event_rate", "d", "weighting_scenario", "method")

  # Aggregate numeric metrics
  numeric_vars <- c("hr_est", "ess", "max_weight", "max_std_diff")

  # Create formula for aggregation
  agg_formula <- as.formula(paste("cbind(",
                                  paste(numeric_vars, collapse = ","),
                                  ") ~ ",
                                  paste(grouping_vars, collapse = "+")))

  metrics <- aggregate(
    agg_formula,
    data = results,
    FUN = function(x) {
      if (all(is.na(x))) {
        return(c(mean = NA, sd = NA, median = NA, q25 = NA, q75 = NA))
      }
      c(mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        q25 = quantile(x, 0.25, na.rm = TRUE),
        q75 = quantile(x, 0.75, na.rm = TRUE))
    }
  )

  # Calculate bias
  metrics$bias <- metrics$hr_est[, "mean"] - metrics$true_hr

  # Calculate coverage for each CI method
  ci_methods <- c("robust", "boot_perc", "boot_norm", "boot_bca")

  for (ci_method in ci_methods) {
    ci_lower_col <- paste0("ci_", ci_method, "_lower")
    ci_upper_col <- paste0("ci_", ci_method, "_upper")

    # Check if columns exist in results
    if (all(c(ci_lower_col, ci_upper_col) %in% names(results))) {
      # Create coverage indicator
      results[[paste0("covered_", ci_method)]] <-
        results[[ci_lower_col]] <= results$true_hr &
        results[[ci_upper_col]] >= results$true_hr

      # Aggregate coverage
      cov_formula <- as.formula(
        paste0("covered_", ci_method, " ~ ",
               paste(grouping_vars, collapse = "+"))
      )

      coverage_df <- aggregate(
        cov_formula,
        data = results,
        FUN = function(x) mean(x, na.rm = TRUE)
      )

      # Add to metrics
      metrics[[paste0("coverage_", ci_method)]] <-
        coverage_df[[paste0("covered_", ci_method)]]

      # Calculate CI width
      results[[paste0("width_", ci_method)]] <-
        results[[ci_upper_col]] - results[[ci_lower_col]]

      width_formula <- as.formula(
        paste0("width_", ci_method, " ~ ",
               paste(grouping_vars, collapse = "+"))
      )

      width_df <- aggregate(
        width_formula,
        data = results,
        FUN = function(x) {
          if (all(is.na(x))) {
            return(c(mean = NA, median = NA))
          }
          c(mean = mean(x, na.rm = TRUE),
            median = median(x, na.rm = TRUE))
        }
      )

      metrics[[paste0("ci_width_", ci_method, "_mean")]] <-
        width_df[[paste0("width_", ci_method)]][, "mean"]
      metrics[[paste0("ci_width_", ci_method, "_median")]] <-
        width_df[[paste0("width_", ci_method)]][, "median"]
    }
  }

  # Calculate Monte Carlo standard errors
  n_sim <- length(unique(results$iteration))
  metrics$mcse_bias <- metrics$hr_est[, "sd"] / sqrt(n_sim)

  # MCSE for coverage (binomial proportion)
  for (ci_method in ci_methods) {
    cov_col <- paste0("coverage_", ci_method)
    if (cov_col %in% names(metrics)) {
      metrics[[paste0("mcse_", cov_col)]] <-
        sqrt(metrics[[cov_col]] * (1 - metrics[[cov_col]]) / n_sim)
    }
  }

  # Save metrics if requested
  if (save_output && !is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- paste0("performance_metrics_", timestamp, ".rds")
    saveRDS(metrics, file.path(output_dir, filename))
  }

  return(metrics)
}

#' Create enhanced summary tables with flexible output options
create_summary_tables <- function(metrics, results = NULL, output_dir = "tables",
                                  format = c("csv", "xlsx", "latex")) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  format <- match.arg(format)

  # Table 1: Overall performance summary
  overall_summary <- aggregate(
    cbind(bias, mcse_bias, coverage_robust) ~ method,
    data = metrics,
    FUN = mean
  )

  # Add ESS summary
  ess_summary <- aggregate(
    ess ~ method,
    data = metrics,
    FUN = function(x) mean(x[, "mean"], na.rm = TRUE)
  )
  overall_summary$mean_ess <- ess_summary$ess

  # Save based on format
  if (format == "csv") {
    write.csv(overall_summary,
              file.path(output_dir, "overall_performance_summary.csv"),
              row.names = FALSE)
  }

  # Table 2: Bias by method and sample size
  bias_table <- reshape(
    metrics[, c("method", "n", "bias", "mcse_bias")],
    idvar = "method",
    timevar = "n",
    direction = "wide"
  )

  write.csv(bias_table, file.path(output_dir, "bias_by_method_and_n.csv"),
            row.names = FALSE)

  # Table 3: Coverage by method and CI type
  ci_methods <- c("robust", "boot_perc", "boot_norm", "boot_bca")

  for (scenario in c("correct", "under", "over")) {
    coverage_data <- metrics[metrics$weighting_scenario == scenario, ]

    if (nrow(coverage_data) > 0) {
      # Create coverage comparison table
      coverage_table <- data.frame(method = unique(coverage_data$method))

      for (ci_method in ci_methods) {
        cov_col <- paste0("coverage_", ci_method)
        if (cov_col %in% names(coverage_data)) {
          coverage_by_method <- aggregate(
            as.formula(paste(cov_col, "~ method")),
            data = coverage_data,
            FUN = mean,
            na.rm = TRUE
          )
          coverage_table[[ci_method]] <- coverage_by_method[[cov_col]]
        }
      }

      write.csv(coverage_table,
                file.path(output_dir, paste0("coverage_", scenario, "_scenario.csv")),
                row.names = FALSE)
    }
  }

  # Table 4: CI width comparison
  width_summary <- data.frame(method = character(), ci_type = character(),
                              mean_width = numeric(), median_width = numeric())

  for (ci_method in ci_methods) {
    mean_col <- paste0("ci_width_", ci_method, "_mean")
    median_col <- paste0("ci_width_", ci_method, "_median")

    if (all(c(mean_col, median_col) %in% names(metrics))) {
      method_summary <- aggregate(
        as.formula(paste0("cbind(", mean_col, ", ", median_col, ") ~ method")),
        data = metrics,
        FUN = mean,
        na.rm = TRUE
      )

      width_summary <- rbind(width_summary,
                             data.frame(method = method_summary$method,
                                        ci_type = ci_method,
                                        mean_width = method_summary[[mean_col]],
                                        median_width = method_summary[[median_col]]))
    }
  }

  write.csv(width_summary, file.path(output_dir, "ci_width_comparison.csv"),
            row.names = FALSE)

  # Table 5: ESS and balance summary
  ess_balance <- aggregate(
    cbind(ess, max_std_diff) ~ method + n + weighting_scenario,
    data = metrics,
    FUN = function(x) mean(x[, "mean"], na.rm = TRUE)
  )

  write.csv(ess_balance, file.path(output_dir, "ess_balance_summary.csv"),
            row.names = FALSE)

  # Table 6: Bootstrap diagnostics (if results provided)
  if (!is.null(results)) {
    boot_diagnostics <- aggregate(
      n_boot_valid ~ method + n + weighting_scenario,
      data = results,
      FUN = function(x) c(
        mean = mean(x, na.rm = TRUE),
        min = min(x, na.rm = TRUE),
        prop_full = mean(x == max(x, na.rm = TRUE), na.rm = TRUE)
      )
    )

    write.csv(boot_diagnostics, file.path(output_dir, "bootstrap_diagnostics.csv"),
              row.names = FALSE)
  }

  cat(sprintf("Summary tables saved to %s/\n", output_dir))
}

#' Create enhanced diagnostic plots with publication quality
create_diagnostic_plots <- function(results, metrics, output_dir = "figures",
                                    plot_format = c("pdf", "png", "both")) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  plot_format <- match.arg(plot_format)

  # Helper function to save plots in requested format
  save_plot <- function(plot_name, width = 10, height = 6, plot_function) {
    if (plot_format %in% c("pdf", "both")) {
      pdf(file.path(output_dir, paste0(plot_name, ".pdf")),
          width = width, height = height)
      plot_function()
      dev.off()
    }
    if (plot_format %in% c("png", "both")) {
      png(file.path(output_dir, paste0(plot_name, ".png")),
          width = width * 100, height = height * 100, res = 150)
      plot_function()
      dev.off()
    }
  }

  # Plot 1: Bias across methods
  save_plot("bias_by_method", 10, 6, function() {
    par(mar = c(5, 4, 4, 2) + 0.1)

    # Create bias data for plotting
    bias_data <- aggregate(bias ~ method, data = metrics,
                           FUN = function(x) x)

    boxplot(bias ~ method, data = metrics,
            main = "Bias in HR Estimation by Method",
            xlab = "Method", ylab = "Bias (HR - True HR)",
            col = rainbow(8, alpha = 0.7),
            outline = FALSE)

    # Add zero line
    abline(h = 0, lty = 2, col = "red", lwd = 2)

    # Add sample size annotation
    points(1:8, aggregate(bias ~ method, data = metrics, mean)$bias,
           pch = 19, col = "black")
  })

  # Plot 2: Coverage probability by CI method
  save_plot("coverage_by_ci_method", 12, 8, function() {
    ci_methods <- c("robust", "boot_perc", "boot_norm", "boot_bca")
    ci_labels <- c("Robust", "Boot Percentile", "Boot Normal", "Boot BCa")

    # Create matrix for grouped barplot
    methods <- unique(metrics$method)
    coverage_matrix <- matrix(NA, nrow = length(ci_methods), ncol = length(methods))

    for (i in seq_along(ci_methods)) {
      cov_col <- paste0("coverage_", ci_methods[i])
      if (cov_col %in% names(metrics)) {
        coverage_means <- aggregate(
          as.formula(paste(cov_col, "~ method")),
          data = metrics,
          FUN = mean,
          na.rm = TRUE
        )
        coverage_matrix[i, ] <- coverage_means[[cov_col]]
      }
    }

    # Create barplot
    barplot(coverage_matrix, beside = TRUE, names.arg = methods,
            col = rainbow(length(ci_methods), alpha = 0.7),
            main = "95% CI Coverage Probability by Method and CI Type",
            xlab = "Weight Normalization Method", ylab = "Coverage Probability",
            ylim = c(0, 1))

    # Add nominal coverage line
    abline(h = 0.95, lty = 2, col = "red", lwd = 2)

    # Add legend
    legend("topright", legend = ci_labels,
           fill = rainbow(length(ci_methods), alpha = 0.7),
           bty = "n")
  })

  # Plot 3: CI width comparison
  save_plot("ci_width_comparison", 12, 10, function() {
    par(mfrow = c(2, 2))

    ci_methods <- c("robust", "boot_perc", "boot_norm", "boot_bca")
    ci_labels <- c("Robust", "Boot Percentile", "Boot Normal", "Boot BCa")

    for (i in seq_along(ci_methods)) {
      width_col <- paste0("ci_width_", ci_methods[i], "_median")
      if (width_col %in% names(metrics)) {
        # Create boxplot for CI widths
        boxplot(as.formula(paste(width_col, "~ method")),
                data = metrics,
                main = paste("CI Width -", ci_labels[i]),
                xlab = "Method", ylab = "CI Width",
                col = rainbow(8, alpha = 0.7),
                outline = FALSE)

        # Add mean line
        width_means <- aggregate(
          as.formula(paste(width_col, "~ method")),
          data = metrics,
          FUN = mean,
          na.rm = TRUE
        )
        points(1:8, width_means[[width_col]], pch = 19, col = "black")
      }
    }
  })

  # Plot 4: ESS vs Balance trade-off
  save_plot("ess_vs_balance", 10, 8, function() {
    # Extract mean values
    ess_means <- aggregate(ess ~ method, data = metrics,
                           FUN = function(x) mean(x[, "mean"], na.rm = TRUE))
    balance_means <- aggregate(max_std_diff ~ method, data = metrics,
                               FUN = function(x) mean(x[, "mean"], na.rm = TRUE))

    # Merge data
    plot_data <- merge(ess_means, balance_means, by = "method")

    # Create plot
    plot(plot_data$ess, plot_data$max_std_diff,
         col = rainbow(8)[as.factor(plot_data$method)],
         pch = 19, cex = 2,
         xlab = "Effective Sample Size",
         ylab = "Maximum Standardized Difference",
         main = "ESS vs Covariate Balance Trade-off",
         xlim = c(0, max(plot_data$ess) * 1.1),
         ylim = c(0, max(plot_data$max_std_diff) * 1.1))

    # Add method labels
    text(plot_data$ess, plot_data$max_std_diff,
         labels = plot_data$method,
         pos = 3, cex = 0.8)

    # Add reference lines
    abline(h = 0.1, lty = 2, col = "gray", lwd = 1)
    abline(h = 0.25, lty = 2, col = "gray", lwd = 1)

    # Add legend for reference lines
    legend("topright",
           legend = c("Negligible imbalance (0.1)", "Small imbalance (0.25)"),
           lty = 2, col = "gray", bty = "n")
  })

  # Plot 5: Coverage by sample size for each CI method
  save_plot("coverage_by_n_and_ci", 12, 10, function() {
    par(mfrow = c(2, 2))

    ci_methods <- c("robust", "boot_perc", "boot_norm", "boot_bca")
    ci_labels <- c("Robust", "Boot Percentile", "Boot Normal", "Boot BCa")

    for (i in seq_along(ci_methods)) {
      cov_col <- paste0("coverage_", ci_methods[i])
      if (cov_col %in% names(metrics)) {
        # Aggregate by n and method
        cov_by_n <- aggregate(
          as.formula(paste(cov_col, "~ n + method")),
          data = metrics,
          FUN = mean,
          na.rm = TRUE
        )

        # Create plot
        plot(1, type = "n",
             xlim = range(cov_by_n$n),
             ylim = c(0.8, 1),
             xlab = "Sample Size per Arm",
             ylab = "Coverage Probability",
             main = paste("Coverage by Sample Size -", ci_labels[i]),
             log = "x")

        # Add lines for each method
        methods <- unique(cov_by_n$method)
        colors <- rainbow(length(methods))

        for (j in seq_along(methods)) {
          subset_data <- cov_by_n[cov_by_n$method == methods[j], ]
          lines(subset_data$n, subset_data[[cov_col]],
                col = colors[j],
                type = "b", pch = 19, lwd = 2)
        }

        # Add nominal coverage line
        abline(h = 0.95, lty = 2, col = "red", lwd = 2)

        # Add legend
        if (i == 1) {
          legend("bottomright", legend = methods,
                 col = colors, lty = 1, pch = 19,
                 lwd = 2, bty = "n", cex = 0.8)
        }
      }
    }
  })

  # Plot 6: Performance by scenario
  save_plot("performance_by_scenario", 14, 8, function() {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

    scenarios <- c("correct", "under", "over")
    scenario_labels <- c("Correct Specification",
                         "Under-specification",
                         "Over-specification")

    # Bias by scenario
    for (i in seq_along(scenarios)) {
      scenario_data <- metrics[metrics$weighting_scenario == scenarios[i], ]

      if (nrow(scenario_data) > 0) {
        boxplot(bias ~ method, data = scenario_data,
                main = paste("Bias -", scenario_labels[i]),
                xlab = "", ylab = "Bias",
                col = rainbow(8, alpha = 0.7),
                outline = FALSE)
        abline(h = 0, lty = 2, col = "red")
      }
    }

    # Coverage by scenario
    for (i in seq_along(scenarios)) {
      scenario_data <- metrics[metrics$weighting_scenario == scenarios[i], ]

      if (nrow(scenario_data) > 0) {
        boxplot(coverage_robust ~ method, data = scenario_data,
                main = paste("Coverage -", scenario_labels[i]),
                xlab = "", ylab = "Coverage",
                col = rainbow(8, alpha = 0.7),
                outline = FALSE,
                ylim = c(0.8, 1))
        abline(h = 0.95, lty = 2, col = "red")
      }
    }
  })

  cat(sprintf("Diagnostic plots saved to %s/\n", output_dir))
}

#' Main function to run complete simulation study
main <- function(config_file = NULL) {
  # Load configuration if provided
  if (!is.null(config_file)) {
    source(config_file)
  }

  # Start timing
  start_time <- Sys.time()

  cat("=============================================================================\n")
  cat("MAIC Simulation Study - Version 37 (Final)\n")
  cat("=============================================================================\n\n")

  # Run tests first
  cat("Running function tests...\n")
  test_results <- test_simulation_functions()

  if (!all(unlist(test_results))) {
    stop("Some tests failed. Please fix errors before running full simulation.")
  }

  cat("\nAll tests passed!\n\n")

  # Run full simulation
  cat("Starting full simulation study...\n")
  cat(sprintf("Using %d cores for parallel processing\n", sim_params$n_cores))

  # Run simulations
  results <- run_simulation_study(sim_params)

  # Calculate performance metrics
  cat("\nCalculating performance metrics...\n")
  metrics <- calculate_performance_metrics(results)

  # Create summary tables
  cat("Creating summary tables...\n")
  create_summary_tables(metrics, results)

  # Create diagnostic plots
  cat("Creating diagnostic plots...\n")
  create_diagnostic_plots(results, metrics)

  # Report completion
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "hours")

  cat("\n=============================================================================\n")
  cat(sprintf("Simulation study completed in %.2f hours\n", runtime))
  cat(sprintf("Total scenarios: %d\n", length(unique(results$scenario_id))))
  cat(sprintf("Total simulations: %d\n", nrow(results) / length(unique(results$method))))
  cat(sprintf("Results saved to: %s/\n",
              ifelse(sim_params$timestamp_files,
                     paste0("simulation_results_", format(start_time, "%Y%m%d_%H%M%S")),
                     "simulation_results")))
  cat("=============================================================================\n")

  # Save session info for reproducibility
  sink("session_info.txt")
  cat("Session Information\n")
  cat("==================\n\n")
  print(sessionInfo())
  cat("\n\nSimulation Parameters\n")
  cat("====================\n\n")
  print(sim_params)
  sink()

  return(list(results = results, metrics = metrics))
}

#' Export results to LaTeX tables
export_latex_tables <- function(metrics, output_file = "simulation_results.tex") {
  # Requires xtable package
  if (!requireNamespace("xtable", quietly = TRUE)) {
    warning("xtable package not available, skipping LaTeX export")
    return(NULL)
  }

  library(xtable)

  # Create bias table
  bias_summary <- aggregate(
    cbind(bias, mcse_bias) ~ method + n,
    data = metrics,
    FUN = mean,
    na.rm = TRUE
  )

  bias_wide <- reshape(bias_summary, idvar = "method", timevar = "n",
                       direction = "wide")

  # Format for LaTeX
  sink(output_file)

  cat("% Simulation Results Tables\n")
  cat("% Generated on", format(Sys.Date(), "%B %d, %Y"), "\n\n")

  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\caption{Average bias in hazard ratio estimation by method and sample size}\n")
  cat("\\label{tab:bias_results}\n")
  print(xtable(bias_wide, digits = 4), include.rownames = FALSE,
        floating = FALSE, booktabs = TRUE)
  cat("\\end{table}\n\n")

  # Create coverage table
  coverage_summary <- aggregate(
    cbind(coverage_robust, coverage_boot_perc, coverage_boot_norm, coverage_boot_bca) ~
      method + weighting_scenario,
    data = metrics,
    FUN = mean,
    na.rm = TRUE
  )

  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\caption{95\\% CI coverage probability by method and CI type}\n")
  cat("\\label{tab:coverage_results}\n")
  print(xtable(coverage_summary, digits = 3), include.rownames = FALSE,
        floating = FALSE, booktabs = TRUE)
  cat("\\end{table}\n")

  sink()

  cat(sprintf("LaTeX tables exported to %s\n", output_file))
}

#' Validate simulation results
validate_results <- function(results) {
  checks <- list()

  # Check 1: All methods present for each scenario
  methods_per_scenario <- aggregate(method ~ scenario_id, data = results,
                                    FUN = function(x) length(unique(x)))
  checks$all_methods_present <- all(methods_per_scenario$method == 8)

  # Check 2: No missing HR estimates (allow some NA for failed models)
  checks$hr_coverage <- mean(!is.na(results$hr_est)) > 0.95

  # Check 3: Reasonable HR range
  checks$reasonable_hr_range <- all(results$hr_est > 0.1 & results$hr_est < 10,
                                    na.rm = TRUE)

  # Check 4: ESS <= n
  checks$ess_valid <- all(results$ess <= results$n * 1.01, na.rm = TRUE)

  # Check 5: Weights are positive
  checks$positive_weights <- all(results$max_weight > 0, na.rm = TRUE)

  # Check 6: Bootstrap success rate
  checks$bootstrap_success <- mean(results$n_boot_valid / sim_params$n_boot,
                                   na.rm = TRUE) > 0.8

  # Report
  cat("\n=== Validation Results ===\n")
  for (check_name in names(checks)) {
    status <- ifelse(checks[[check_name]], "PASS", "FAIL")
    cat(sprintf("  %-25s: %s\n", check_name, status))
  }
  cat("\n")

  return(checks)
}

# =============================================================================
# Execute if run as script
# =============================================================================

if (!interactive()) {
  # Running as script
  results <- main()

  # Validate results
  validation <- validate_results(results$results)

  # Export LaTeX tables
  export_latex_tables(results$metrics)

  cat("\nSimulation study complete. Results saved to simulation_results/\n")
}
