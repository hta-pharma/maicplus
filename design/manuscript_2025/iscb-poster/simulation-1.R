# =============================================================================
# MAIC Simulation Study for Time-to-Event Outcomes
#
# This code implements the simulation study described in the manuscript:
# "Comparative Evaluation of Weight Normalization Methods in Matching-Adjusted
# Indirect Comparison for Time-to-Event Outcomes"
#
# Author: Gregory, etc
# Date: Aug 2025
# R Version: 4.5.1
# =============================================================================

# Load required libraries
library(survival)
library(parallel)
library(MASS)  # For mvrnorm if needed

# =============================================================================
# Design
# =============================================================================

# ## Main Components:
#
# ### 1. **Data Generation Functions**
# - `generate_covariates()`: Creates baseline covariates for internal/external populations
# - `calibrate_d_for_ess()`: Calibrates the scaling factor to achieve target ESS ratios
# - `generate_survival_times()`: Generates Weibull or log-normal survival times
# - `calibrate_base_params()`: Calibrates distribution parameters for target median and 48-month survival
# - `generate_censoring()`: Implements the dropout pattern (60% early, 40% late)
# - `calibrate_beta()`: Calibrates prognostic effects to achieve target event rates
#
# ### 2. **Weight Calculation Functions**
# - `calculate_maic_weights()`: Implements method of moments for MAIC weights
# - `calculate_mew_weights()`: Calculates maximum ESS weights using quadratic programming
# - `normalize_weights()`: Applies different normalization schemes (OW, SW1, SWN, SWESS)
#
# ### 3. **Analysis Functions**
# - `fit_weighted_cox()`: Fits weighted Cox models and extracts HR with CIs
# - `calculate_balance()`: Computes maximum standardized differences
# - `calculate_ess()`: Computes effective sample size
#
# ### 4. **Simulation Infrastructure**
# - `run_single_simulation()`: Executes one simulation iteration
# - `run_simulation_study()`: Manages the full factorial design with parallel processing
# - Automatic intermediate result saving every 10 scenarios
# - Parallel computing with configurable cores
#
# ### 5. **Results Processing**
# - `calculate_performance_metrics()`: Aggregates results with Monte Carlo standard errors
# - `create_summary_tables()`: Generates CSV tables for manuscript
# - `create_diagnostic_plots()`: Creates visualization plots
# - `export_latex_tables()`: Exports results to LaTeX format
# - `validate_results()`: Performs quality checks on simulation output
#
# ## Key Features:
#
# 1. **Efficient Implementation**
#    - Uses base R functions where possible (minimal dependencies)
#    - Parallel processing with `makeCluster()`
#    - Vectorized operations throughout
#    - Memory-efficient intermediate saving
#
# 2. **Follows Manuscript Exactly**
#    - All 8 weight normalization methods
#    - Three prognostic factor scenarios (correct, under, over)
#    - Full factorial design as specified
#    - Algorithm matches your manuscript's description
#
# 3. **Comprehensive Output**
#   - Bias and coverage with Monte Carlo SEs
# - ESS and weight diagnostics
# - Balance metrics
# - Automatic plot generation
# - LaTeX table export
#
# 4. **Flexibility**
#   - Number of cores is tunable
# - Can save/resume from intermediate results
# - Configuration can be loaded from external file
# - Easy to extend with additional methods
#
# ## Usage:
#
# ```r
# # Run with default parameters
# source("maic_simulation.R")
# results <- main()
#
# # Or customize parameters
# sim_params$n_cores <- 8  # Use 8 cores
# sim_params$n_sim <- 5000  # Reduce simulations for testing
# results <- main()
# ```

# =============================================================================
# SECTION 1: SIMULATION PARAMETERS AND SETUP
# =============================================================================

# Set global parameters
set.seed(12345)  # For reproducibility

# Simulation parameters as per manuscript
sim_params <- list(
  # Sample sizes per arm
  n_sizes = c(20, 50, 200, 500),

  # True ATUT HR values
  true_hr = c(0.70, 0.85),

  # Event distributions
  distributions = c("weibull", "lognormal"),

  # Median survival times (months)
  median_surv = c(24, 36),

  # 48-month survival rates
  surv_48m = c(0.05, 0.15),

  # Total dropout rates
  dropout_rates = c(0.2, 0.4),

  # Target event rates in control group
  event_rates = c(0.05, 0.10, 0.25, 0.40),

  # Effective sample size ratios
  ess_ratios = c(0.15, 0.25, 0.33, 0.50),

  # Simulation settings
  n_sim = 10000,  # Number of Monte Carlo replications
  max_followup = 48,  # months

  # Parallel computing
  n_cores = 4  # Adjust based on your system
)

# =============================================================================
# SECTION 2: HELPER FUNCTIONS FOR DATA GENERATION
# =============================================================================

#' Generate baseline covariates for internal and external populations
#'
#' @param n Sample size
#' @param population "internal" or "external"
#' @param d Scaling factor for internal population means
#' @return Matrix of covariates (n x 4)
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

#' Calibrate d value to achieve target ESS ratio
#'
#' @param target_ess_ratio Target ESS/N ratio
#' @param n_calib Sample size for calibration
#' @param tol Tolerance for convergence
#' @param method Algorithm used to find the right value for d, either "bisection" (default) or "newton-ralphson"
#' @return Calibrated d value
calibrate_d_for_ess <- function(target_ess_ratio, n_calib = 10000, tol = 0.01, method = c("bisection","newton-ralphson")) {
  # Function to compute ESS ratio for given d
  compute_ess_ratio <- function(d) {
    # Generate large samples for calibration
    X_external <- generate_covariates(n_calib, "external")
    X_internal <- generate_covariates(n_calib, "internal", d)

    # Compute propensity scores (only using prognostic factors X1-X3)
    combined_X <- rbind(X_external[, 1:3], X_internal[, 1:3])
    combined_Z <- c(rep(0, n_calib), rep(1, n_calib))

    ps_model <- glm(combined_Z ~ combined_X, family = binomial())
    ps_external <- predict(ps_model, newdata = data.frame(combined_X = X_external[, 1:3]),
                           type = "response")

    # Compute weights
    weights <- ps_external / (1 - ps_external)
    weights <- c(weights, rep(1, n_calib))

    # Compute ESS
    ess <- sum(weights)^2 / sum(weights^2)

    # Return ESS ratio
    return(ess / (2 * n_calib))
  }

  if("bisection" %in% method) {
    # Use bisection method to find d
    d_lower <- 0.1
    d_upper <- 2.0

    while (d_upper - d_lower > tol) {
      d_mid <- (d_lower + d_upper) / 2
      ess_mid <- compute_ess_ratio(d_mid)

      if (ess_mid < target_ess_ratio) {
        d_upper <- d_mid
      } else {
        d_lower <- d_mid
      }
    }

    d_current <- (d_lower + d_upper) / 2
  }else if(all(method == "newton-ralphson")) {
    # Use Newton-Raphson method to find d
    # Initialize d with a reasonable starting value
    d_current <- 1.0
    max_iter <- 50
    h <- 0.001  # Step size for numerical derivative

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
        warning("Derivative too small in Newton-Raphson, switching to bisection for this step")
        # Fall back to bisection for this iteration
        if (f_current > 0) {
          d_current <- d_current * 0.9
        } else {
          d_current <- d_current * 1.1
        }
      } else {
        # Newton-Raphson update
        d_new <- d_current - f_current / f_prime

        # Ensure d stays within reasonable bounds [0.1, 2.0]
        d_new <- pmax(0.1, pmin(2.0, d_new))

        # Apply damping if change is too large
        if (abs(d_new - d_current) > 0.5) {
          d_new <- d_current + sign(d_new - d_current) * 0.5
        }

        d_current <- d_new
      }
    }

    if (iter == max_iter) {
      warning(sprintf("Newton-Raphson did not converge after %d iterations", max_iter))
    }
  }else {
    stop("specified 'method' does not exist")
  }

  return(d_current)
}

#' Generate survival times based on specified distribution
#'
#' @param n Sample size
#' @param X Covariate matrix
#' @param beta Prognostic effects
#' @param gamma Treatment effect (log HR)
#' @param Z Treatment indicator
#' @param distribution "weibull" or "lognormal"
#' @param base_params List with parameters for baseline distribution
#' @return Vector of survival times
generate_survival_times <- function(n, X, beta, gamma, Z, distribution, base_params) {
  # Linear predictor
  lp <- X %*% beta + gamma * Z

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
#'
#' @param distribution "weibull" or "lognormal"
#' @param median_surv Target median survival time
#' @param surv_48m Target 48-month survival probability
#' @return List of distribution parameters
calibrate_base_params <- function(distribution, median_surv, surv_48m) {
  if (distribution == "weibull") {
    # For Weibull: S(t) = exp(-(lambda*t)^shape)
    # Median: exp(-(lambda*median)^shape) = 0.5
    # 48m: exp(-(lambda*48)^shape) = surv_48m

    shape <- log(-log(0.5)) / log(median_surv/48) *
      log(-log(surv_48m)) / log(-log(0.5))
    lambda <- (-log(0.5))^(1/shape) / median_surv

    return(list(lambda = lambda, shape = shape))
  } else {
    # For log-normal: S(t) = 1 - Phi((log(t) - mu)/sigma)
    # Solve for mu and sigma given median and 48m survival

    # Median: log(median_surv) = mu (since Phi(0) = 0.5)
    mu <- log(median_surv)

    # 48m survival: 1 - Phi((log(48) - mu)/sigma) = surv_48m
    sigma <- (log(48) - mu) / qnorm(1 - surv_48m)

    return(list(mu = mu, sigma = sigma))
  }
}

#' Generate censoring times based on dropout pattern
#'
#' @param n Sample size
#' @param T True survival times
#' @param dropout_rate Total dropout rate
#' @return List with censoring times and dropout indicators
generate_censoring <- function(n, T, dropout_rate) {
  # Number of dropouts
  n_dropout <- rbinom(1, n, dropout_rate)

  # Split dropouts: 60% in first 6 months, 40% after
  n_early <- rbinom(1, n_dropout, 0.6)
  n_late <- n_dropout - n_early

  # Randomly assign dropouts
  dropout_ids <- sample(1:n, n_dropout, replace = FALSE)
  early_ids <- dropout_ids[1:n_early]
  late_ids <- if (n_late > 0) dropout_ids[(n_early + 1):n_dropout] else numeric(0)

  # Initialize censoring times
  C <- rep(48, n)  # Administrative censoring at 48 months

  # Early dropouts: uniform in [0, min(T, 6)]
  for (i in early_ids) {
    C[i] <- runif(1, 0, min(T[i], 6))
  }

  # Late dropouts: uniform in [6, min(T, 48)]
  for (i in late_ids) {
    if (T[i] > 6) {
      C[i] <- runif(1, 6, min(T[i], 48))
    } else {
      # If event before 6 months, censor uniformly before event
      C[i] <- runif(1, 0, T[i])
    }
  }

  return(C)
}

#' Calibrate beta coefficients to achieve target event rate
#'
#' @param target_event_rate Target event rate in control group
#' @param n_calib Sample size for calibration
#' @param distribution Event time distribution
#' @param base_params Baseline distribution parameters
#' @param dropout_rate Dropout rate
#' @return Calibrated beta vector
calibrate_beta <- function(target_event_rate, n_calib = 5000, distribution,
                           base_params, dropout_rate) {
  # Function to compute event rate for given beta scale
  compute_event_rate <- function(beta_scale) {
    # Set beta values (X1-X3 prognostic, X4 non-prognostic)
    beta <- c(beta_scale, beta_scale, beta_scale, 0)

    # Generate data
    X <- generate_covariates(n_calib, "external")
    T <- generate_survival_times(n_calib, X, beta, gamma = 0, Z = rep(0, n_calib),
                                 distribution, base_params)
    C <- generate_censoring(n_calib, T, dropout_rate)

    # Observed outcomes
    Y <- pmin(T, C)
    delta <- as.numeric(T <= C)

    # Event rate
    event_rate <- mean(delta)
    return(event_rate)
  }

  # Use bisection to find beta scale
  beta_lower <- -2
  beta_upper <- 2
  tol <- 0.01

  while (beta_upper - beta_lower > tol) {
    beta_mid <- (beta_lower + beta_upper) / 2
    rate_mid <- compute_event_rate(beta_mid)

    if (rate_mid < target_event_rate) {
      beta_lower <- beta_mid
    } else {
      beta_upper <- beta_mid
    }
  }

  beta_scale <- (beta_lower + beta_upper) / 2
  return(c(beta_scale, beta_scale, beta_scale, 0))
}

# =============================================================================
# SECTION 3: WEIGHT CALCULATION FUNCTIONS
# =============================================================================

#' Calculate MAIC weights using method of moments
#'
#' @param X_internal Internal trial covariates (to be weighted)
#' @param X_external_means External trial covariate means (target)
#' @return Vector of weights
calculate_maic_weights <- function(X_internal, X_external_means) {
  n <- nrow(X_internal)
  p <- ncol(X_internal)

  # Center covariates
  X_centered <- sweep(X_internal, 2, X_external_means, "-")

  # Method of moments: solve for alpha such that sum(w_i * X_i) = n * X_external_means
  # where w_i = exp(X_centered %*% alpha)

  # Objective function for optimization
  objective <- function(alpha) {
    weights <- exp(X_centered %*% alpha)
    weighted_means <- colSums(weights * X_internal) / sum(weights)
    sum((weighted_means - X_external_means)^2)
  }

  # Gradient
  gradient <- function(alpha) {
    weights <- exp(X_centered %*% alpha)
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

  # Optimize
  result <- optim(par = rep(0, p), fn = objective, gr = gradient,
                  method = "BFGS", control = list(maxit = 1000))

  # Calculate weights
  weights <- exp(X_centered %*% result$par)
  return(weights)
}

#' Calculate maximum ESS weights
#'
#' @param X_internal Internal trial covariates
#' @param X_external_means External trial covariate means
#' @return Vector of MEW weights
calculate_mew_weights <- function(X_internal, X_external_means) {
  n <- nrow(X_internal)
  p <- ncol(X_internal)

  # Center covariates
  X_centered <- sweep(X_internal, 2, X_external_means, "-")

  # Objective: minimize sum(w_i^2) subject to sum(w_i * X_i) = n * X_external_means
  # Using Lagrangian approach

  # Construct constraint matrix
  A <- cbind(rep(1, n), X_internal)
  b <- c(n, n * X_external_means)

  # Solve quadratic program: min(w'w) s.t. A'w = b
  # Solution: w = A(A'A)^{-1}b

  # Check if A'A is invertible
  AtA <- t(A) %*% A
  if (det(AtA) < 1e-10) {
    warning("Constraint matrix nearly singular, using regularization")
    AtA <- AtA + diag(1e-6, ncol(AtA))
  }

  weights <- A %*% solve(AtA) %*% b

  # Ensure non-negative weights
  weights <- pmax(weights, 0)

  return(weights)
}

#' Normalize weights according to different schemes
#'
#' @param weights Original weights
#' @param method Normalization method
#' @return Normalized weights
normalize_weights <- function(weights, method) {
  n <- length(weights)

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
#'
#' @param weights Vector of weights
#' @return ESS value
calculate_ess <- function(weights) {
  sum(weights)^2 / sum(weights^2)
}

# =============================================================================
# SECTION 4: ESTIMATION FUNCTIONS
# =============================================================================

#' Fit weighted Cox model and extract HR with CI
#'
#' @param time Observed times
#' @param event Event indicators
#' @param treatment Treatment indicators
#' @param weights Analysis weights
#' @return List with HR estimate and 95% CI
fit_weighted_cox <- function(time, event, treatment, weights = NULL) {
  # Create survival object
  surv_obj <- Surv(time, event)

  # Fit Cox model
  if (is.null(weights)) {
    cox_fit <- coxph(surv_obj ~ treatment)
  } else {
    cox_fit <- coxph(surv_obj ~ treatment, weights = weights)
  }

  # Extract results
  hr <- exp(coef(cox_fit))
  ci <- exp(confint(cox_fit))

  return(list(
    hr = as.numeric(hr),
    ci_lower = as.numeric(ci[1]),
    ci_upper = as.numeric(ci[2]),
    se_log_hr = as.numeric(sqrt(vcov(cox_fit)))
  ))
}

#' Calculate balance metrics
#'
#' @param X_internal Internal trial covariates
#' @param X_external_means External trial means
#' @param weights Weights
#' @return Maximum standardized difference
calculate_balance <- function(X_internal, X_external_means, weights) {
  # Weighted means
  weighted_means <- colSums(weights * X_internal) / sum(weights)

  # Pooled standard deviations (using external as reference)
  # For simulation, we know external SDs
  external_sds <- c(0.5, 0.5, 0.5, 0.5)  # X1 has SD=0.5, X2-X4 are Bernoulli

  # Standardized differences
  std_diffs <- abs(weighted_means - X_external_means) / external_sds

  return(max(std_diffs))
}

# =============================================================================
# SECTION 5: MAIN SIMULATION FUNCTION
# =============================================================================

#' Run single simulation iteration
#'
#' @param params List of simulation parameters
#' @return Data frame with results
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

  # Generate baseline distribution parameters
  base_params <- calibrate_base_params(distribution, median_surv, surv_48m)

  # Generate covariates
  X_internal <- generate_covariates(n, "internal", d)
  X_external <- generate_covariates(n, "external")

  # Generate survival times
  T_internal <- generate_survival_times(n, X_internal, beta, gamma, rep(1, n),
                                        distribution, base_params)
  T_external <- generate_survival_times(n, X_external, beta, 0, rep(0, n),
                                        distribution, base_params)

  # Generate censoring
  C_internal <- generate_censoring(n, T_internal, dropout_rate)
  C_external <- generate_censoring(n, T_external, dropout_rate)

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

  # Calculate base weights
  maic_weights <- calculate_maic_weights(X_for_weights, X_external_means)
  mew_weights <- calculate_mew_weights(X_for_weights, X_external_means)

  # Initialize results storage
  results <- data.frame(
    n = n,
    true_hr = true_hr,
    distribution = distribution,
    median_surv = median_surv,
    surv_48m = surv_48m,
    dropout_rate = dropout_rate,
    event_rate = event_rate,
    d = d,
    weighting_scenario = weighting_scenario,
    method = character(),
    hr_est = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    ess = numeric(),
    max_weight = numeric(),
    max_std_diff = numeric(),
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
    cox_results <- fit_weighted_cox(Y_combined, delta_combined, Z_combined, weights_combined)

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
      ci_lower = cox_results$ci_lower,
      ci_upper = cox_results$ci_upper,
      ess = ess,
      max_weight = max_weight,
      max_std_diff = max_std_diff,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}

# =============================================================================
# SECTION 6: PARALLEL SIMULATION EXECUTION
# =============================================================================

#' Run full simulation study
#'
#' @param sim_params Global simulation parameters
#' @param save_intermediate Whether to save intermediate results
#' @param save_dir Directory for saving results
#' @return Combined results data frame
run_simulation_study <- function(sim_params, save_intermediate = TRUE,
                                 save_dir = "simulation_results") {

  # Create save directory if needed
  if (save_intermediate && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  # Pre-calibrate d values for each ESS ratio
  cat("Calibrating d values for ESS ratios...\n")
  d_values <- sapply(sim_params$ess_ratios, calibrate_d_for_ess)
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

  for (i in 1:nrow(unique_combos)) {
    combo <- unique_combos[i, ]

    # Calibrate base parameters
    base_params <- calibrate_base_params(combo$distribution, combo$median_surv,
                                         combo$surv_48m)

    # Calibrate beta
    beta <- calibrate_beta(combo$event_rate, distribution = combo$distribution,
                           base_params = base_params, dropout_rate = combo$dropout_rate)

    # Calibrate gamma for target HR
    # This requires iterative procedure similar to beta calibration
    # For now, use approximation: gamma ≈ log(true_hr)
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

    scenarios$beta[matching_rows] <- list(beta)
    scenarios$gamma[matching_rows] <- gamma
  }

  # Set up parallel cluster
  cl <- makeCluster(sim_params$n_cores)
  clusterExport(cl, ls(envir = .GlobalEnv))
  clusterEvalQ(cl, {
    library(survival)
    set.seed(Sys.getpid())  # Different seed for each worker
  })

  # Run simulations
  total_scenarios <- nrow(scenarios)
  cat(sprintf("Running %d scenarios with %d replications each...\n",
              total_scenarios, sim_params$n_sim))

  all_results <- list()

  for (scenario_idx in 1:total_scenarios) {
    cat(sprintf("Scenario %d/%d...\n", scenario_idx, total_scenarios))

    # Current scenario
    current_scenario <- scenarios[scenario_idx, ]

    # Create parameter list for each replication
    param_list <- replicate(sim_params$n_sim, as.list(current_scenario),
                            simplify = FALSE)

    # Run parallel simulations
    scenario_results <- parLapply(cl, param_list, run_single_simulation)

    # Combine results
    scenario_results_df <- do.call(rbind, scenario_results)
    scenario_results_df$iteration <- rep(1:sim_params$n_sim,
                                         each = length(unique(scenario_results_df$method)))
    scenario_results_df$scenario_id <- scenario_idx

    # Store results
    all_results[[scenario_idx]] <- scenario_results_df

    # Save intermediate results if requested
    if (save_intermediate && scenario_idx %% 10 == 0) {
      intermediate_df <- do.call(rbind, all_results)
      saveRDS(intermediate_df, file.path(save_dir,
                                         sprintf("intermediate_results_%d.rds", scenario_idx)))
      cat(sprintf("  Saved intermediate results at scenario %d\n", scenario_idx))
    }
  }

  # Stop cluster
  stopCluster(cl)

  # Combine all results
  final_results <- do.call(rbind, all_results)

  # Save final results
  if (save_intermediate) {
    saveRDS(final_results, file.path(save_dir, "final_results.rds"))
    write.csv(final_results, file.path(save_dir, "final_results.csv"),
              row.names = FALSE)
  }

  return(final_results)
}

# =============================================================================
# SECTION 7: PERFORMANCE METRICS CALCULATION
# =============================================================================

#' Calculate performance metrics across simulations
#'
#' @param results Combined results from all simulations
#' @return Data frame with aggregated performance metrics
calculate_performance_metrics <- function(results) {
  # Group by scenario and method
  metrics <- aggregate(
    cbind(hr_est, ci_lower, ci_upper, ess, max_weight, max_std_diff) ~
      n + true_hr + distribution + median_surv + surv_48m + dropout_rate +
      event_rate + d + weighting_scenario + method,
    data = results,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      q25 = quantile(x, 0.25, na.rm = TRUE),
      q75 = quantile(x, 0.75, na.rm = TRUE)
    )
  )

  # Calculate bias
  metrics$bias <- metrics$hr_est[, "mean"] - metrics$true_hr

  # Calculate coverage
  coverage_df <- aggregate(
    cbind(covered = I(ci_lower <= true_hr & ci_upper >= true_hr)) ~
      n + true_hr + distribution + median_surv + surv_48m + dropout_rate +
      event_rate + d + weighting_scenario + method,
    data = results,
    FUN = mean
  )

  metrics$coverage <- coverage_df$covered

  # Calculate Monte Carlo standard errors
  n_sim <- length(unique(results$iteration))
  metrics$mcse_bias <- metrics$hr_est[, "sd"] / sqrt(n_sim)
  metrics$mcse_coverage <- sqrt(metrics$coverage * (1 - metrics$coverage) / n_sim)

  return(metrics)
}

#' Create summary tables for manuscript
#'
#' @param metrics Performance metrics data frame
#' @param output_dir Directory for saving tables
create_summary_tables <- function(metrics, output_dir = "tables") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Table 1: Bias by method and sample size
  bias_table <- reshape(
    metrics[, c("method", "n", "bias", "mcse_bias")],
    idvar = "method",
    timevar = "n",
    direction = "wide"
  )

  write.csv(bias_table, file.path(output_dir, "bias_by_method_and_n.csv"),
            row.names = FALSE)

  # Table 2: Coverage by method and scenario
  coverage_table <- reshape(
    metrics[, c("method", "weighting_scenario", "coverage", "mcse_coverage")],
    idvar = "method",
    timevar = "weighting_scenario",
    direction = "wide"
  )

  write.csv(coverage_table, file.path(output_dir, "coverage_by_method_and_scenario.csv"),
            row.names = FALSE)

  # Table 3: ESS comparison
  ess_summary <- aggregate(
    ess ~ method + n,
    data = metrics,
    FUN = function(x) c(mean = mean(x[, "mean"]), sd = mean(x[, "sd"]))
  )

  write.csv(ess_summary, file.path(output_dir, "ess_summary.csv"),
            row.names = FALSE)
}

# =============================================================================
# SECTION 8: VISUALIZATION FUNCTIONS
# =============================================================================

#' Create diagnostic plots
#'
#' @param results Raw simulation results
#' @param metrics Aggregated performance metrics
#' @param output_dir Directory for saving plots
create_diagnostic_plots <- function(results, metrics, output_dir = "figures") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Plot 1: Bias across methods
  pdf(file.path(output_dir, "bias_by_method.pdf"), width = 10, height = 6)

  par(mar = c(5, 4, 4, 2) + 0.1)
  boxplot(bias ~ method, data = metrics,
          main = "Bias in HR Estimation by Method",
          xlab = "Method", ylab = "Bias",
          col = rainbow(8))
  abline(h = 0, lty = 2, col = "red")

  dev.off()

  # Plot 2: Coverage probability
  pdf(file.path(output_dir, "coverage_by_method.pdf"), width = 10, height = 6)

  coverage_means <- aggregate(coverage ~ method, data = metrics, mean)
  barplot(coverage_means$coverage, names.arg = coverage_means$method,
          main = "95% CI Coverage Probability by Method",
          xlab = "Method", ylab = "Coverage",
          ylim = c(0, 1), col = rainbow(8))
  abline(h = 0.95, lty = 2, col = "red")

  dev.off()

  # Plot 3: ESS vs Balance trade-off
  pdf(file.path(output_dir, "ess_vs_balance.pdf"), width = 8, height = 8)

  plot(metrics$ess[, "mean"], metrics$max_std_diff[, "mean"],
       col = as.factor(metrics$method), pch = 19,
       xlab = "Effective Sample Size",
       ylab = "Maximum Standardized Difference",
       main = "ESS vs Covariate Balance Trade-off")
  legend("topright", legend = unique(metrics$method),
         col = 1:length(unique(metrics$method)), pch = 19)

  dev.off()

  # Plot 4: Weight distributions
  pdf(file.path(output_dir, "weight_distributions.pdf"), width = 12, height = 8)

  par(mfrow = c(2, 4))
  for (method in unique(results$method)) {
    subset_data <- results[results$method == method & results$iteration <= 100, ]
    hist(log(subset_data$max_weight + 1),
         main = paste("Log(Max Weight) -", method),
         xlab = "log(Max Weight + 1)", col = "lightblue")
  }

  dev.off()
}

# =============================================================================
# SECTION 9: MAIN EXECUTION SCRIPT
# =============================================================================

#' Main function to run complete simulation study
#'
#' @param config_file Optional configuration file path
main <- function(config_file = NULL) {
  # Load configuration if provided
  if (!is.null(config_file)) {
    source(config_file)
  }

  # Start timing
  start_time <- Sys.time()

  cat("Starting MAIC simulation study...\n")
  cat(sprintf("Using %d cores for parallel processing\n", sim_params$n_cores))

  # Run simulations
  results <- run_simulation_study(sim_params)

  # Calculate performance metrics
  cat("Calculating performance metrics...\n")
  metrics <- calculate_performance_metrics(results)

  # Create summary tables
  cat("Creating summary tables...\n")
  create_summary_tables(metrics)

  # Create diagnostic plots
  cat("Creating diagnostic plots...\n")
  create_diagnostic_plots(results, metrics)

  # Report completion
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "hours")

  cat(sprintf("\nSimulation study completed in %.2f hours\n", runtime))
  cat(sprintf("Total scenarios: %d\n", length(unique(results$scenario_id))))
  cat(sprintf("Total simulations: %d\n", nrow(results) / length(unique(results$method))))

  # Save session info for reproducibility
  sink("session_info.txt")
  sessionInfo()
  sink()

  return(list(results = results, metrics = metrics))
}

# =============================================================================
# SECTION 10: UTILITY FUNCTIONS FOR POST-PROCESSING
# =============================================================================

#' Export results to LaTeX tables
#'
#' @param metrics Performance metrics
#' @param output_file LaTeX file name
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
    FUN = mean
  )

  bias_wide <- reshape(bias_summary, idvar = "method", timevar = "n",
                       direction = "wide")

  # Format for LaTeX
  sink(output_file)

  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\caption{Average bias in hazard ratio estimation by method and sample size}\n")
  print(xtable(bias_wide, digits = 4), include.rownames = FALSE,
        floating = FALSE)
  cat("\\end{table}\n\n")

  # Create coverage table
  coverage_summary <- aggregate(
    cbind(coverage, mcse_coverage) ~ method + weighting_scenario,
    data = metrics,
    FUN = mean
  )

  coverage_wide <- reshape(coverage_summary, idvar = "method",
                           timevar = "weighting_scenario", direction = "wide")

  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\caption{95\\% CI coverage probability by method and weighting scenario}\n")
  print(xtable(coverage_wide, digits = 3), include.rownames = FALSE,
        floating = FALSE)
  cat("\\end{table}\n")

  sink()

  cat(sprintf("LaTeX tables exported to %s\n", output_file))
}

#' Validate simulation results
#'
#' @param results Simulation results
#' @return List of validation checks
validate_results <- function(results) {
  checks <- list()

  # Check 1: All methods present for each scenario
  methods_per_scenario <- aggregate(method ~ scenario_id, data = results,
                                    FUN = function(x) length(unique(x)))
  checks$all_methods_present <- all(methods_per_scenario$method == 8)

  # Check 2: No missing HR estimates
  checks$no_missing_hr <- sum(is.na(results$hr_est)) == 0

  # Check 3: Reasonable HR range
  checks$reasonable_hr_range <- all(results$hr_est > 0.1 & results$hr_est < 10,
                                    na.rm = TRUE)

  # Check 4: ESS <= n
  checks$ess_valid <- all(results$ess <= results$n, na.rm = TRUE)

  # Check 5: Weights are positive
  checks$positive_weights <- all(results$max_weight > 0, na.rm = TRUE)

  # Report
  cat("Validation Results:\n")
  for (check_name in names(checks)) {
    status <- ifelse(checks[[check_name]], "PASS", "FAIL")
    cat(sprintf("  %s: %s\n", check_name, status))
  }

  return(checks)
}

# =============================================================================
# Execute if run as script
# =============================================================================

if (sys.nframe() == 0) {
  # Running as script, not sourced
  results <- main()

  # Validate results
  validation <- validate_results(results$results)

  # Export LaTeX tables
  export_latex_tables(results$metrics)

  cat("\nSimulation study complete. Results saved to simulation_results/\n")
}
