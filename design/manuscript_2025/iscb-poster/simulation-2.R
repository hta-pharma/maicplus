# =============================================================================
# MAIC Simulation Study for Time-to-Event Outcomes
#
# This code implements the simulation study described in the manuscript:
# "Comparative Evaluation of Weight Normalization Methods in Matching-Adjusted
# Indirect Comparison for Time-to-Event Outcomes"
#
# Author: Gregory Chen, etc
# Date: Aug 10, 2025
# R Version: 4.5.1
# =============================================================================

# =============================================================================
# HIGH-LEVEL DESIGN SUMMARY
# =============================================================================
#
# This simulation study evaluates different weight normalization methods for
# Matching-Adjusted Indirect Comparison (MAIC) in time-to-event analyses.
#
# MAIN COMPONENTS:
#
# 1. DATA GENERATION (Section 2)
#    - generate_covariates(): Creates baseline characteristics
#    - generate_survival_times(): Simulates event times (Weibull/Log-normal)
#    - generate_censoring(): Implements realistic dropout patterns
#    - calibrate_d_for_ess(): Calibrates population differences
#    - calibrate_beta(): Calibrates prognostic effects
#    - calibrate_base_params(): Sets baseline hazard parameters
#
# 2. WEIGHT CALCULATION (Section 3)
#    - calculate_maic_weights(): Standard MAIC weights via method of moments
#    - calculate_mew_weights(): Maximum ESS weights via optimization
#    - normalize_weights(): 8 different normalization schemes
#
# 3. ANALYSIS METHODS (Section 4)
#    - fit_weighted_cox(): Cox regression with multiple CI methods
#    - calculate_balance(): Covariate balance assessment
#    - calculate_ess(): Effective sample size computation
#
# 4. SIMULATION FRAMEWORK (Sections 5-6)
#    - run_single_simulation(): Executes one complete simulation
#    - run_simulation_study(): Manages factorial design with parallelization
#
# 5. RESULTS PROCESSING (Sections 7-10)
#    - calculate_performance_metrics(): Aggregates results
#    - create_summary_tables(): Generates publication-ready tables
#    - create_diagnostic_plots(): Creates visualization figures
#    - export_latex_tables(): Formats results for manuscripts
#
# EXECUTION WORKFLOW:
#
# # Step 1: Load the script
# source("maic_simulation.R")
#
# # Step 2: Configure simulation parameters (optional)
# sim_params$n_cores <- 8        # Use 8 CPU cores
# sim_params$n_sim <- 1000       # Reduce for testing (default 10,000)
# sim_params$n_boot <- 200       # Reduce bootstrap iterations (default 500)
#
# # Step 3: Run the full simulation study
# results <- main()
#
# # Step 4: Access results
# # results$results contains raw simulation data
# # results$metrics contains aggregated performance metrics
#
# # Step 5: Generate additional outputs
# export_latex_tables(results$metrics, "manuscript_tables.tex")
#
# # For testing with a subset:
# test_params <- sim_params
# test_params$n_sizes <- c(50, 200)
# test_params$n_sim <- 100
# test_params$true_hr <- 0.7
# test_results <- run_simulation_study(test_params, save_intermediate = FALSE)
#
# =============================================================================

# Load required libraries
library(survival)
library(parallel)
library(MASS)  # For mvrnorm if needed

# =============================================================================
# SECTION 1: SIMULATION PARAMETERS AND SETUP
# =============================================================================

# Set global parameters
set.seed(12345)  # For reproducibility

# Simulation parameters as per manuscript Table 1
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
  n_cores = 4  # Adjust based on your system
)

# =============================================================================
# SECTION 2: HELPER FUNCTIONS FOR DATA GENERATION
# =============================================================================

#' Generate baseline covariates for internal and external populations
#'
#' @description Creates four baseline covariates (X1-X4) with controlled
#' differences between internal and external populations. X1-X3 are prognostic
#' factors affecting survival, while X4 is a non-prognostic nuisance variable.
#'
#' @param n Integer. Sample size
#' @param population Character. Either "internal" or "external"
#' @param d Numeric. Scaling factor for internal population means (1 = no difference).
#'   Controls the degree of covariate imbalance between populations.
#'
#' @return Matrix of dimensions n x 4 containing:
#'   \item{X1}{Continuous covariate ~ N(0.5*d, 0.5^2)}
#'   \item{X2}{Binary covariate ~ Bernoulli(min(0.5*d, 1))}
#'   \item{X3}{Binary covariate ~ Bernoulli(min(0.5*d, 1))}
#'   \item{X4}{Binary covariate ~ Bernoulli(min(0.5*d, 1)) (non-prognostic)}
#'
#' @details The external population has fixed parameters with X1 ~ N(0.5, 0.5^2)
#' and X2-X4 ~ Bernoulli(0.5). The internal population parameters are scaled by d,
#' creating systematic differences that affect the effective sample size when
#' matching populations via weighting.
#'
#' @examples
#' # Generate balanced populations (d = 1)
#' X_ext <- generate_covariates(100, "external", d = 1)
#' X_int <- generate_covariates(100, "internal", d = 1)
#'
#' # Generate imbalanced populations (d = 1.5)
#' X_ext <- generate_covariates(100, "external", d = 1)
#' X_int <- generate_covariates(100, "internal", d = 1.5)
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
#' @description Uses optimization to find the scaling factor d that produces
#' a desired effective sample size ratio when reweighting the external population
#' to match the internal population distribution.
#'
#' @param target_ess_ratio Numeric in (0,1). Target ESS/N ratio where N is total
#'   sample size (internal + external populations)
#' @param n_calib Integer. Sample size for calibration (default 10,000).
#'   Larger values give more accurate calibration but take longer.
#' @param tol Numeric. Tolerance for convergence (default 0.01)
#' @param method Character. Algorithm: "bisection" (default) or "newton-ralphson"
#'
#' @return Numeric. Calibrated d value that achieves target ESS ratio
#'
#' @details This function solves for d such that when external population is
#' reweighted to match internal population (with mean vector d*mu_0), the
#' resulting ESS/N equals the target. Uses propensity score weighting with
#' weights = e/(1-e) where e is propensity score from logistic regression.
#'
#' The relationship between d and ESS is U-shaped with maximum at d=1 (no
#' difference). This implementation searches for solutions where d>1 to ensure
#' consistency between methods.
#'
#' @references
#' Signorovitch et al. (2010). Comparative effectiveness without head-to-head
#' trials. Pharmacoeconomics 28(10):935-945.
#'
#' @examples
#' # Find d for 25% ESS ratio (high imbalance)
#' d_high <- calibrate_d_for_ess(0.25)
#'
#' # Find d for 50% ESS ratio (moderate imbalance)
#' d_mod <- calibrate_d_for_ess(0.50)
calibrate_d_for_ess <- function(target_ess_ratio, n_calib = 10000, tol = 0.01, method = c("bisection","newton-ralphson")) {
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

    ps_model <- glm(Z ~ X1 + X2 + X3, data = glm_data, family = binomial())

    # Create prediction data frame
    pred_data <- data.frame(
      X1 = X_external[, 1],
      X2 = X_external[, 2],
      X3 = X_external[, 3]
    )

    ps_external <- predict(ps_model, newdata = pred_data, type = "response")

    # Compute weights with numerical stability
    # Avoid extreme propensity scores
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
    # Since ESS is highest when d = 1 and decreases as d moves away from 1,
    # we search in [1, 2] for the solution where d > 1
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
    # Start from d slightly above 1 to ensure we find the d > 1 solution
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
        warning("Derivative too small in Newton-Raphson, switching to bisection for this step")
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

    if (iter == max_iter) {
      warning(sprintf("Newton-Raphson did not converge after %d iterations", max_iter))
    }
  } else {
    stop("specified 'method' does not exist")
  }

  return(d_current)
}

#' Generate survival times based on specified distribution
#'
#' @description Generates event times from either Weibull or log-normal
#' distributions with covariates affecting the hazard via Cox regression model.
#'
#' @param n Integer. Sample size
#' @param X Matrix. Covariate matrix (n x p)
#' @param beta Numeric vector. Regression coefficients for covariates
#' @param gamma Numeric. Log hazard ratio for treatment effect
#' @param Z Numeric vector. Treatment indicators (0/1)
#' @param distribution Character. Either "weibull" or "lognormal"
#' @param base_params List. Parameters for baseline distribution:
#'   For Weibull: list(lambda = scale, shape = shape)
#'   For log-normal: list(mu = mean, sigma = sd)
#'
#' @return Numeric vector of survival times
#'
#' @details Implements the proportional hazards model:
#'   h(t|X,Z) = h0(t) * exp(X*beta + gamma*Z)
#'
#' For Weibull: Uses inverse CDF method with
#'   T = (-log(U) / (lambda * exp(X*beta + gamma*Z)))^(1/shape)
#'
#' For log-normal: Generates from
#'   log(T) ~ N(mu - X*beta - gamma*Z, sigma^2)
#'
#' @references
#' Bender et al. (2005). Generating survival times to simulate Cox proportional
#' hazards models. Statistics in Medicine 24(11):1713-1723.
#'
#' @examples
#' # Weibull with median 24 months
#' base_params <- list(lambda = 0.03, shape = 1.5)
#' X <- matrix(rnorm(100*3), 100, 3)
#' T <- generate_survival_times(100, X, beta = c(0.5, -0.3, 0.2),
#'                             gamma = log(0.7), Z = rep(0:1, each = 50),
#'                             "weibull", base_params)
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
#' @description Calculates distribution parameters to achieve specified median
#' survival time and survival probability at a target time point.
#'
#' @param distribution Character. "weibull" or "lognormal"
#' @param median_surv Numeric. Target median survival time
#' @param surv_target Numeric. Target survival probability at target_time
#' @param target_time Numeric. Time point for survival probability (default 48)
#'
#' @return List of distribution parameters:
#'   For Weibull: list(lambda = scale, shape = shape)
#'   For log-normal: list(mu = location, sigma = scale)
#'
#' @details Solves the system of equations:
#'   S(median_surv) = 0.5
#'   S(target_time) = surv_target
#'
#' For Weibull: S(t) = exp(-(lambda*t)^shape)
#' For log-normal: S(t) = 1 - Phi((log(t) - mu)/sigma)
#'
#' @examples
#' # Weibull with median 24 months and 15% survival at 48 months
#' params_wb <- calibrate_base_params("weibull", 24, 0.15, 48)
#'
#' # Log-normal with same targets
#' params_ln <- calibrate_base_params("lognormal", 24, 0.15, 48)
calibrate_base_params <- function(distribution, median_surv, surv_target, target_time = 48) {
  if (distribution == "weibull") {
    # For Weibull: S(t) = exp(-(lambda*t)^shape)
    # Median: exp(-(lambda*median)^shape) = 0.5
    # Target time: exp(-(lambda*target_time)^shape) = surv_target

    shape <- log(-log(0.5)) / log(median_surv/target_time) *
      log(-log(surv_target)) / log(-log(0.5))
    lambda <- (-log(0.5))^(1/shape) / median_surv

    return(list(lambda = lambda, shape = shape))
  } else {
    # For log-normal: S(t) = 1 - Phi((log(t) - mu)/sigma)
    # Solve for mu and sigma given median and target_time survival

    # Median: log(median_surv) = mu (since Phi(0) = 0.5)
    mu <- log(median_surv)

    # Target time survival: 1 - Phi((log(target_time) - mu)/sigma) = surv_target
    sigma <- (log(target_time) - mu) / qnorm(1 - surv_target)

    return(list(mu = mu, sigma = sigma))
  }
}

#' Generate censoring times based on dropout pattern
#'
#' @description Implements realistic dropout patterns as specified in the
#' simulation design. Ensures exactly the specified dropout rate with 60%
#' occurring in first 6 months.
#'
#' @param n Integer. Sample size
#' @param T Numeric vector. True survival times
#' @param dropout_rate Numeric. Proportion who dropout (excluding admin censoring)
#' @param max_followup Numeric. Maximum follow-up time (default 48 months)
#'
#' @return Numeric vector of censoring times
#'
#' @details Implements the dropout mechanism where:
#' - Exactly dropout_rate proportion will dropout before end of study
#' - 60% of dropouts occur uniformly in [0, 6] months
#' - 40% of dropouts occur uniformly in [6, max_followup] months
#' - Dropout times are generated conditional on true event times to ensure
#'   dropout occurs before the event would have been observed
#' - All others are administratively censored at max_followup
#'
#' @examples
#' # 20% dropout rate
#' T <- rexp(100, rate = 0.02)
#' C <- generate_censoring(100, T, 0.2, 48)
#'
#' # Check dropout occurred before events
#' all(C[C < 48] < T[C < 48])  # Should be TRUE
generate_censoring <- function(n, T, dropout_rate, max_followup = 48) {
  # Number of dropouts
  n_dropout <- rbinom(1, n, dropout_rate)

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

#' Calibrate beta coefficients to achieve target event rate
#'
#' @description Finds the prognostic effect size (beta) that produces a desired
#' event rate in the control group, accounting for censoring.
#'
#' @param target_event_rate Numeric. Target event rate in control group
#' @param n_calib Integer. Sample size for calibration (default 5000)
#' @param distribution Character. Event time distribution
#' @param base_params List. Baseline distribution parameters
#' @param dropout_rate Numeric. Dropout rate
#' @param max_followup Numeric. Maximum follow-up time (default 48)
#' @param tol Numeric. Tolerance for convergence (default 0.01)
#' @param method Character. "bisection" (default) or "newton-ralphson"
#'
#' @return Numeric vector. Beta coefficients where first 3 are equal (prognostic)
#'   and fourth is 0 (non-prognostic)
#'
#' @details Uses optimization to find beta such that the observed event rate
#' (accounting for dropout and administrative censoring) equals the target.
#' Assumes equal effects for prognostic factors X1-X3 and no effect for X4.
#'
#' @examples
#' base_params <- calibrate_base_params("weibull", 24, 0.15)
#' beta <- calibrate_beta(0.25, distribution = "weibull",
#'                       base_params = base_params,
#'                       dropout_rate = 0.2)
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

    # Verify monotonic increasing relationship
    if (rate_lower >= rate_upper) {
      warning("Non-monotonic relationship detected between beta and event rate")
    }

    # Check if target is achievable
    if (target_event_rate < rate_lower || target_event_rate > rate_upper) {
      warning(sprintf("Target event rate %.3f outside achievable range [%.3f, %.3f]",
                      target_event_rate, rate_lower, rate_upper))
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
      warning(sprintf("Target event rate %.3f outside achievable range [%.3f, %.3f]",
                      target_event_rate, rate_at_minus2, rate_at_plus2))
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
        warning("Derivative too small in Newton-Raphson, switching to bisection for this step")
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

    if (iter == max_iter) {
      warning(sprintf("Newton-Raphson did not converge after %d iterations", max_iter))
    }

    beta_scale <- beta_current
  } else {
    stop("specified 'method' does not exist")
  }

  return(c(beta_scale, beta_scale, beta_scale, 0))
}

# =============================================================================
# SECTION 3: WEIGHT CALCULATION FUNCTIONS
# =============================================================================

#' Calculate MAIC weights using method of moments
#'
#' @description Computes balancing weights for the internal (treated) population
#' to match the external (control) population's covariate distribution using
#' the method of moments approach.
#'
#' @param X_internal Matrix. Internal trial covariates to be weighted (n x p)
#' @param X_external_means Numeric vector. External trial covariate means (length p)
#'
#' @return Numeric vector of weights (length n)
#'
#' @details Solves for weights w_i = exp(X_centered_i * alpha) such that:
#'   sum(w_i * X_i) / sum(w_i) = X_external_means
#'
#' where X_centered = X_internal - X_external_means
#'
#' Uses BFGS optimization with analytical gradient. The objective function
#' minimizes squared differences between weighted means and target means.
#'
#' @references
#' Signorovitch et al. (2010). Comparative effectiveness without head-to-head
#' trials. Pharmacoeconomics 28(10):935-945.
#'
#' Phillippo et al. (2018). Methods for population-adjusted indirect comparisons
#' in health technology appraisal. Medical Decision Making 38(2):200-211.
#'
#' @examples
#' X_int <- matrix(rnorm(100*3, mean = 0.5), 100, 3)
#' X_ext_means <- c(0, 0, 0)
#' weights <- calculate_maic_weights(X_int, X_ext_means)
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
#' @description Computes weights that maximize effective sample size (ESS)
#' subject to exact balance constraints. This provides optimal efficiency
#' while achieving covariate balance.
#'
#' @param X_internal Matrix. Internal trial covariates (n x p)
#' @param X_external_means Numeric vector. External trial covariate means
#'
#' @return Numeric vector of MEW weights
#'
#' @details Solves the optimization problem:
#'   minimize sum(w_i^2)
#'   subject to: sum(w_i) = n
#'               sum(w_i * X_i) = n * X_external_means
#'
#' This is equivalent to maximizing ESS = (sum w_i)^2 / sum(w_i^2)
#'
#' Uses closed-form solution: w = A(A'A)^{-1}b where A = [1, X] and
#' b = [n, n*X_external_means]
#'
#' @references
#' Zubizarreta (2015). Stable weights that balance covariates for estimation
#' with incomplete outcome data. JASA 110(511):910-922.
#'
#' @examples
#' X_int <- matrix(rnorm(100*3, mean = 0.5), 100, 3)
#' X_ext_means <- c(0, 0, 0)
#' mew_weights <- calculate_mew_weights(X_int, X_ext_means)
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
#' @description Applies various normalization schemes to raw weights, affecting
#' the interpretation and performance of weighted estimators.
#'
#' @param weights Numeric vector. Original weights
#' @param method Character. Normalization method:
#'   - "OW": Original weights (unscaled)
#'   - "SW1": Sum to 1
#'   - "SWN": Sum to N (sample size)
#'   - "SWESS": Sum to ESS (effective sample size)
#'
#' @return Numeric vector of normalized weights
#'
#' @details Different normalizations affect:
#' - Variance estimation (robust vs model-based)
#' - Interpretation (average vs total)
#' - Small sample properties
#'
#' ESS = (sum w_i)^2 / sum(w_i^2) measures the effective number of
#' independent observations after weighting.
#'
#' @examples
#' w <- c(0.5, 1.0, 1.5, 2.0)
#' w_norm1 <- normalize_weights(w, "SW1")    # sum to 1
#' w_normN <- normalize_weights(w, "SWN")    # sum to N
#' w_normESS <- normalize_weights(w, "SWESS") # sum to ESS
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
#' @description Computes the effective sample size (ESS) after weighting,
#' which quantifies the loss of precision due to variable weights.
#'
#' @param weights Numeric vector of weights
#'
#' @return Numeric. ESS value
#'
#' @details ESS = (sum w_i)^2 / sum(w_i^2)
#'
#' Interpretation:
#' - ESS = n: all weights equal (no loss of precision)
#' - ESS < n: unequal weights reduce precision
#' - ESS/n: proportion of effective information retained
#'
#' @references
#' Kish (1965). Survey Sampling. New York: Wiley.
#'
#' @examples
#' # Equal weights: ESS = n
#' calculate_ess(rep(1, 100))  # Returns 100
#'
#' # Unequal weights: ESS < n
#' calculate_ess(c(rep(0.5, 50), rep(1.5, 50)))  # Returns < 100
calculate_ess <- function(weights) {
  sum(weights)^2 / sum(weights^2)
}

# =============================================================================
# SECTION 4: ESTIMATION FUNCTIONS
# =============================================================================

#' Fit weighted Cox model and extract HR with multiple CI methods
#'
#' @description Fits a weighted Cox proportional hazards model and computes
#' the hazard ratio with multiple confidence interval methods for comparison.
#'
#' @param time Numeric vector. Observed times
#' @param event Logical/numeric vector. Event indicators (1 = event, 0 = censored)
#' @param treatment Numeric vector. Treatment indicators (1 = treated, 0 = control)
#' @param weights Numeric vector. Analysis weights (NULL for unweighted)
#' @param n_boot Integer. Number of bootstrap iterations (default 500)
#'
#' @return List containing:
#'   \item{hr}{Point estimate of hazard ratio}
#'   \item{log_hr}{Log hazard ratio}
#'   \item{se_log_hr}{Robust standard error of log HR}
#'   \item{ci_robust_lower/upper}{Robust variance CI (Lin-Wei)}
#'   \item{ci_boot_perc_lower/upper}{Bootstrap percentile CI}
#'   \item{ci_boot_norm_lower/upper}{Bootstrap normal approximation CI}
#'   \item{ci_boot_bca_lower/upper}{BCa bootstrap CI}
#'   \item{n_boot_valid}{Number of successful bootstrap iterations}
#'   \item{boot_log_hr_dist}{Bootstrap distribution of log HR}
#'
#' @details Implements four CI methods:
#' 1. Robust: Uses Lin-Wei sandwich variance (always computed)
#' 2. Bootstrap percentile: 2.5th and 97.5th percentiles of bootstrap HRs
#' 3. Bootstrap normal: Point estimate ± 1.96*SE from bootstrap
#' 4. BCa: Bias-corrected and accelerated bootstrap with jackknife
#'
#' Bootstrap uses stratified sampling by treatment group to ensure balance.
#'
#' @references
#' Lin & Wei (1989). The robust inference for the Cox proportional hazards
#' model. JASA 84(408):1074-1078.
#'
#' Efron & Tibshirani (1993). An Introduction to the Bootstrap. Chapman & Hall.
#'
#' @examples
#' # Generate example data
#' n <- 200
#' time <- rexp(n, 0.1)
#' event <- rbinom(n, 1, 0.7)
#' treatment <- rep(0:1, each = n/2)
#' weights <- runif(n, 0.5, 2)
#'
#' # Fit model with all CI methods
#' results <- fit_weighted_cox(time, event, treatment, weights, n_boot = 200)
fit_weighted_cox <- function(time, event, treatment, weights = NULL, n_boot = 500) {
  # Create survival object
  surv_obj <- Surv(time, event)

  # Fit Cox model with robust variance if weights are provided
  if (is.null(weights)) {
    cox_fit <- coxph(surv_obj ~ treatment)
    robust_se <- sqrt(vcov(cox_fit)[1,1])
  } else {
    # Use robust=TRUE for sandwich variance estimator
    cox_fit <- coxph(surv_obj ~ treatment, weights = weights, robust = TRUE)
    robust_se <- sqrt(cox_fit$var[1,1])  # Robust variance
  }

  # Extract point estimate
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
      boot_idx_0 <- sample(trt_0_idx, n_0, replace = TRUE)
      boot_idx_1 <- sample(trt_1_idx, n_1, replace = TRUE)
      boot_idx <- c(boot_idx_0, boot_idx_1)

      # Bootstrap data
      boot_time <- time[boot_idx]
      boot_event <- event[boot_idx]
      boot_treatment <- treatment[boot_idx]
      boot_weights <- if (!is.null(weights)) weights[boot_idx] else NULL

      # Fit model on bootstrap sample
      tryCatch({
        boot_surv <- Surv(boot_time, boot_event)
        if (is.null(boot_weights)) {
          boot_cox <- coxph(boot_surv ~ boot_treatment)
        } else {
          boot_cox <- coxph(boot_surv ~ boot_treatment, weights = boot_weights)
        }
        boot_log_hr[b] <- as.numeric(coef(boot_cox))
        boot_valid[b] <- TRUE
      }, error = function(e) {
        boot_valid[b] <- FALSE
      })
    }

    # Remove invalid bootstrap samples
    boot_log_hr_valid <- boot_log_hr[boot_valid]
    n_valid <- sum(boot_valid)

    if (n_valid >= 50) {  # Need reasonable number for CIs
      # Method 2: Bootstrap percentile CI
      boot_hr <- exp(boot_log_hr_valid)
      ci_boot_perc_lower <- quantile(boot_hr, 0.025, na.rm = TRUE)
      ci_boot_perc_upper <- quantile(boot_hr, 0.975, na.rm = TRUE)

      # Method 3: Bootstrap normal CI (using SE from bootstrap)
      boot_se <- sd(boot_log_hr_valid, na.rm = TRUE)
      ci_boot_norm_lower <- exp(log_hr - 1.96 * boot_se)
      ci_boot_norm_upper <- exp(log_hr + 1.96 * boot_se)

      # Method 4: BCa (Bias-corrected and accelerated) bootstrap CI
      # Calculate bias correction factor z0
      p_less <- mean(boot_log_hr_valid < log_hr)
      z0 <- qnorm(p_less)

      # Calculate acceleration factor using jackknife
      jack_log_hr <- numeric(n_total)
      jack_valid <- logical(n_total)

      # Jackknife to estimate acceleration
      for (i in 1:n_total) {
        jack_idx <- seq_len(n_total)[-i]
        tryCatch({
          jack_surv <- Surv(time[jack_idx], event[jack_idx])
          if (is.null(weights)) {
            jack_cox <- coxph(jack_surv ~ treatment[jack_idx])
          } else {
            jack_cox <- coxph(jack_surv ~ treatment[jack_idx],
                              weights = weights[jack_idx])
          }
          jack_log_hr[i] <- as.numeric(coef(jack_cox))
          jack_valid[i] <- TRUE
        }, error = function(e) {
          jack_valid[i] <- FALSE
        })
      }

      if (sum(jack_valid) >= n_total * 0.8) {  # Need most jackknife samples
        jack_log_hr_valid <- jack_log_hr[jack_valid]
        jack_mean <- mean(jack_log_hr_valid)
        jack_diff <- jack_mean - jack_log_hr_valid

        # Acceleration factor
        a <- sum(jack_diff^3) / (6 * (sum(jack_diff^2))^(3/2))

        # BCa adjusted percentiles
        alpha_lower <- 0.025
        alpha_upper <- 0.975
        z_lower <- qnorm(alpha_lower)
        z_upper <- qnorm(alpha_upper)

        # Adjusted percentiles
        p_lower <- pnorm(z0 + (z0 + z_lower) / (1 - a * (z0 + z_lower)))
        p_upper <- pnorm(z0 + (z0 + z_upper) / (1 - a * (z0 + z_upper)))

        # Ensure percentiles are in [0,1]
        p_lower <- max(0.001, min(0.999, p_lower))
        p_upper <- max(0.001, min(0.999, p_upper))

        # BCa CI
        ci_boot_bca_lower <- quantile(boot_hr, p_lower, na.rm = TRUE)
        ci_boot_bca_upper <- quantile(boot_hr, p_upper, na.rm = TRUE)
      }

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
    n_boot_valid = if(n_boot > 0) sum(boot_valid) else NA,
    boot_log_hr_dist = boot_log_hr_dist
  ))
}

#' Calculate balance metrics
#'
#' @description Computes standardized differences to assess covariate balance
#' after weighting. Returns the maximum absolute standardized difference.
#'
#' @param X_internal Matrix. Internal trial covariates (n x p)
#' @param X_external_means Numeric vector. External trial means (target)
#' @param weights Numeric vector. Weights for internal population
#'
#' @return Numeric. Maximum absolute standardized difference across covariates
#'
#' @details Standardized difference for covariate j:
#'   d_j = |weighted_mean_j - target_mean_j| / sd_j
#'
#' where sd_j is the standard deviation in the external population.
#'
#' Common thresholds:
#' - |d| < 0.1: Negligible imbalance
#' - |d| < 0.25: Small imbalance
#' - |d| > 0.25: Substantial imbalance
#'
#' @references
#' Austin (2009). Balance diagnostics for comparing the distribution of
#' baseline covariates between treatment groups. Statistics in Medicine
#' 28(25):3083-3107.
#'
#' @examples
#' X <- matrix(rnorm(100*3), 100, 3)
#' target_means <- c(0, 0, 0)
#' weights <- runif(100, 0.5, 2)
#' max_diff <- calculate_balance(X, target_means, weights)
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
#' @description Executes one complete simulation iteration including data
#' generation, weight calculation, model fitting, and performance evaluation
#' for all weight normalization methods.
#'
#' @param params List containing all simulation parameters:
#'   - n: Sample size per arm
#'   - true_hr: True ATUT hazard ratio
#'   - distribution: "weibull" or "lognormal"
#'   - median_surv: Median survival time
#'   - surv_48m: 48-month survival probability
#'   - dropout_rate: Dropout proportion
#'   - event_rate: Target event rate
#'   - d: Population difference scaling factor
#'   - beta: Prognostic effect coefficients
#'   - gamma: Log hazard ratio
#'   - weighting_scenario: "correct", "under", or "over"
#'   - max_followup: Maximum follow-up time
#'
#' @return Data frame with one row per weight normalization method containing:
#'   - All input parameters
#'   - HR estimate and all CI bounds
#'   - ESS, maximum weight, balance metrics
#'   - Bootstrap diagnostics
#'
#' @details Workflow:
#' 1. Generate covariates for both populations
#' 2. Generate survival times with treatment effect
#' 3. Apply censoring mechanism
#' 4. Calculate MAIC and MEW weights
#' 5. Apply 8 normalization schemes
#' 6. Fit weighted Cox models with multiple CIs
#' 7. Assess balance and ESS
#'
#' @examples
#' # Set up parameters for one scenario
#' params <- list(
#'   n = 100, true_hr = 0.7, distribution = "weibull",
#'   median_surv = 24, surv_48m = 0.15, dropout_rate = 0.2,
#'   event_rate = 0.25, d = 1.5, beta = c(0.5, 0.5, 0.5, 0),
#'   gamma = log(0.7), weighting_scenario = "correct",
#'   max_followup = 48
#' )
#' results <- run_single_simulation(params)
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
  max_followup <- params$max_followup  # Add this parameter

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
      # Store all CI types
      ci_robust_lower = cox_results$ci_robust_lower,
      ci_robust_upper = cox_results$ci_robust_upper,
      ci_boot_perc_lower = cox_results$ci_boot_perc_lower,
      ci_boot_perc_upper = cox_results$ci_boot_perc_upper,
      ci_boot_norm_lower = cox_results$ci_boot_norm_lower,
      ci_boot_norm_upper = cox_results$ci_boot_norm_upper,
      ci_boot_bca_lower = cox_results$ci_boot_bca_lower,
      ci_boot_bca_upper = cox_results$ci_boot_bca_upper,
      # Other metrics
      ess = ess,
      max_weight = max_weight,
      max_std_diff = max_std_diff,
      n_boot_valid = cox_results$n_boot_valid,
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
#' @description Executes the complete factorial simulation design using parallel
#' processing. Manages scenario generation, parameter calibration, and result
#' aggregation with optional intermediate saving.
#'
#' @param sim_params List. Global simulation parameters including:
#'   - n_sizes: Sample sizes to evaluate
#'   - true_hr: True hazard ratios
#'   - distributions: Event time distributions
#'   - n_sim: Number of replications per scenario
#'   - n_cores: Number of CPU cores for parallel processing
#'   - n_boot: Bootstrap iterations per simulation
#' @param save_intermediate Logical. Save results every 10 scenarios (default TRUE)
#' @param save_dir Character. Directory for saving results (default "simulation_results")
#'
#' @return Data frame containing all simulation results with columns for:
#'   - Scenario parameters
#'   - Method identifiers
#'   - Performance metrics
#'   - Iteration and scenario IDs
#'
#' @details The function:
#' 1. Creates full factorial design of all parameter combinations
#' 2. Pre-calibrates d values for each ESS ratio
#' 3. Calibrates beta and gamma for each unique combination
#' 4. Distributes simulations across CPU cores
#' 5. Saves intermediate results to prevent data loss
#' 6. Combines all results into final dataset
#'
#' Total scenarios = product of all parameter levels × 3 weighting scenarios
#' Total simulations = scenarios × n_sim × 8 methods
#'
#' @examples
#' # Run reduced simulation for testing
#' test_params <- sim_params
#' test_params$n_sizes <- c(50, 200)
#' test_params$n_sim <- 100
#' test_params$n_cores <- 2
#'
#' results <- run_simulation_study(test_params, save_intermediate = FALSE)
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
                                         combo$surv_48m, target_time = sim_params$max_followup)

    # Calibrate beta
    beta <- calibrate_beta(combo$event_rate, distribution = combo$distribution,
                           base_params = base_params, dropout_rate = combo$dropout_rate,
                           max_followup = sim_params$max_followup)

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
    current_scenario$max_followup <- sim_params$max_followup  # Add max_followup

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
#' @description Aggregates raw simulation results to compute performance metrics
#' including bias, coverage, CI width, and Monte Carlo standard errors for each
#' scenario and method combination.
#'
#' @param results Data frame. Combined results from all simulations
#'
#' @return Data frame with aggregated metrics containing:
#'   - Scenario identifiers (n, true_hr, distribution, etc.)
#'   - Method identifier
#'   - Bias and MCSE of bias
#'   - Coverage and MCSE for each CI method
#'   - CI width statistics for each CI method
#'   - ESS and balance summaries
#'
#' @details Computes:
#' - Bias: mean(HR_est) - true_HR
#' - Coverage: proportion of CIs containing true HR
#' - CI width: upper - lower bounds
#' - Monte Carlo SE: SD/sqrt(n_sim) for continuous, sqrt(p(1-p)/n_sim) for proportions
#'
#' @references
#' Morris et al. (2019). Using simulation studies to evaluate statistical
#' methods. Statistics in Medicine 38(11):2074-2102.
#'
#' @examples
#' # Assuming results from simulation
#' metrics <- calculate_performance_metrics(simulation_results)
#'
#' # Extract coverage for robust CI
#' coverage_robust <- metrics$coverage_robust
#'
#' # Compare CI widths
#' width_comparison <- metrics[, grep("ci_width", names(metrics))]
calculate_performance_metrics <- function(results) {
  # Group by scenario and method
  grouping_vars <- c("n", "true_hr", "distribution", "median_surv", "surv_48m",
                     "dropout_rate", "event_rate", "d", "weighting_scenario", "method")

  # Aggregate numeric metrics
  numeric_vars <- c("hr_est", "ess", "max_weight", "max_std_diff")

  metrics <- aggregate(
    as.formula(paste("cbind(", paste(numeric_vars, collapse = ","), ") ~ ",
                     paste(grouping_vars, collapse = "+"))),
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

  # Calculate coverage for each CI method
  ci_methods <- c("robust", "boot_perc", "boot_norm", "boot_bca")

  for (ci_method in ci_methods) {
    ci_lower_col <- paste0("ci_", ci_method, "_lower")
    ci_upper_col <- paste0("ci_", ci_method, "_upper")

    # Check if columns exist in results
    if (all(c(ci_lower_col, ci_upper_col) %in% names(results))) {
      coverage_formula <- as.formula(
        paste0("cbind(covered = I(", ci_lower_col, " <= true_hr & ",
               ci_upper_col, " >= true_hr)) ~ ",
               paste(grouping_vars, collapse = "+"))
      )

      coverage_df <- aggregate(
        coverage_formula,
        data = results,
        FUN = function(x) mean(x, na.rm = TRUE)
      )

      # Add to metrics
      metrics[[paste0("coverage_", ci_method)]] <- coverage_df$covered

      # Calculate CI width
      width_formula <- as.formula(
        paste0("cbind(width = ", ci_upper_col, " - ", ci_lower_col, ") ~ ",
               paste(grouping_vars, collapse = "+"))
      )

      width_df <- aggregate(
        width_formula,
        data = results,
        FUN = function(x) c(
          mean = mean(x, na.rm = TRUE),
          median = median(x, na.rm = TRUE)
        )
      )

      metrics[[paste0("ci_width_", ci_method, "_mean")]] <- width_df$width[, "mean"]
      metrics[[paste0("ci_width_", ci_method, "_median")]] <- width_df$width[, "median"]
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

  return(metrics)
}

#' Create summary tables for manuscript
#'
#' @description Generates publication-ready tables summarizing simulation results,
#' organized by different performance aspects and saved as CSV files.
#'
#' @param metrics Data frame. Aggregated performance metrics
#' @param output_dir Character. Directory for saving tables (default "tables")
#'
#' @details Creates the following tables:
#' 1. Bias by method and sample size
#' 2. Coverage by CI method and weighting scenario
#' 3. CI width comparison across methods
#' 4. ESS summary statistics
#' 5. Bootstrap convergence diagnostics
#'
#' All tables are saved as CSV files for easy import into manuscripts.
#'
#' @examples
#' metrics <- calculate_performance_metrics(results)
#' create_summary_tables(metrics, "manuscript_tables")
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

  # Table 2: Coverage by method and CI type
  ci_methods <- c("robust", "boot_perc", "boot_norm", "boot_bca")

  for (scenario in c("correct", "under", "over")) {
    coverage_data <- metrics[metrics$weighting_scenario == scenario, ]

    # Create coverage comparison table
    coverage_table <- data.frame(method = unique(coverage_data$method))

    for (ci_method in ci_methods) {
      cov_col <- paste0("coverage_", ci_method)
      if (cov_col %in% names(coverage_data)) {
        coverage_by_method <- aggregate(
          as.formula(paste(cov_col, "~ method")),
          data = coverage_data,
          FUN = mean
        )
        coverage_table[[ci_method]] <- coverage_by_method[[cov_col]]
      }
    }

    write.csv(coverage_table,
              file.path(output_dir, paste0("coverage_", scenario, "_scenario.csv")),
              row.names = FALSE)
  }

  # Table 3: CI width comparison
  width_comparison <- data.frame(
    method = unique(metrics$method),
    scenario = rep(unique(metrics$weighting_scenario),
                   each = length(unique(metrics$method)))
  )

  for (ci_method in ci_methods) {
    width_col <- paste0("ci_width_", ci_method, "_median")
    if (width_col %in% names(metrics)) {
      width_data <- aggregate(
        as.formula(paste(width_col, "~ method + weighting_scenario")),
        data = metrics,
        FUN = mean
      )
      width_comparison[[paste0(ci_method, "_width")]] <- width_data[[width_col]]
    }
  }

  write.csv(width_comparison, file.path(output_dir, "ci_width_comparison.csv"),
            row.names = FALSE)

  # Table 4: ESS comparison
  ess_summary <- aggregate(
    ess ~ method + n,
    data = metrics,
    FUN = function(x) c(mean = mean(x[, "mean"]), sd = mean(x[, "sd"]))
  )

  write.csv(ess_summary, file.path(output_dir, "ess_summary.csv"),
            row.names = FALSE)

  # Table 5: Bootstrap diagnostics
  boot_diagnostics <- aggregate(
    n_boot_valid ~ method + n + weighting_scenario,
    data = results,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      prop_full = mean(x == 500, na.rm = TRUE)
    )
  )

  write.csv(boot_diagnostics, file.path(output_dir, "bootstrap_diagnostics.csv"),
            row.names = FALSE)
}

# =============================================================================
# SECTION 8: VISUALIZATION FUNCTIONS
# =============================================================================

#' Create diagnostic plots
#'
#' @description Generates comprehensive visualization of simulation results
#' including bias, coverage, CI width, and diagnostic plots.
#'
#' @param results Data frame. Raw simulation results
#' @param metrics Data frame. Aggregated performance metrics
#' @param output_dir Character. Directory for saving plots (default "figures")
#'
#' @details Creates the following plots:
#' 1. Bias boxplots by method
#' 2. Coverage probability by CI method (grouped barplot)
#' 3. CI width comparison (panel plot)
#' 4. ESS vs balance trade-off (scatter plot)
#' 5. Coverage by sample size for each CI method
#' 6. Bootstrap convergence diagnostics
#'
#' All plots are saved as PDF files with publication-quality formatting.
#'
#' @examples
#' results <- run_simulation_study(sim_params)
#' metrics <- calculate_performance_metrics(results)
#' create_diagnostic_plots(results, metrics, "simulation_figures")
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

  # Plot 2: Coverage probability by CI method
  pdf(file.path(output_dir, "coverage_by_ci_method.pdf"), width = 12, height = 8)

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
        FUN = mean
      )
      coverage_matrix[i, ] <- coverage_means[[cov_col]]
    }
  }

  barplot(coverage_matrix, beside = TRUE, names.arg = methods,
          col = rainbow(length(ci_methods)),
          main = "95% CI Coverage Probability by Method and CI Type",
          xlab = "Weight Normalization Method", ylab = "Coverage",
          ylim = c(0, 1))
  abline(h = 0.95, lty = 2, col = "red")
  legend("topright", legend = ci_labels, fill = rainbow(length(ci_methods)))

  dev.off()

  # Plot 3: CI width comparison
  pdf(file.path(output_dir, "ci_width_comparison.pdf"), width = 12, height = 8)

  par(mfrow = c(2, 2))
  for (i in seq_along(ci_methods)) {
    width_col <- paste0("ci_width_", ci_methods[i], "_median")
    if (width_col %in% names(metrics)) {
      boxplot(as.formula(paste(width_col, "~ method")), data = metrics,
              main = paste("CI Width -", ci_labels[i]),
              xlab = "Method", ylab = "CI Width",
              col = rainbow(8))
    }
  }

  dev.off()

  # Plot 4: ESS vs Balance trade-off
  pdf(file.path(output_dir, "ess_vs_balance.pdf"), width = 8, height = 8)

  plot(metrics$ess[, "mean"], metrics$max_std_diff[, "mean"],
       col = as.factor(metrics$method), pch = 19,
       xlab = "Effective Sample Size",
       ylab = "Maximum Standardized Difference",
       main = "ESS vs Covariate Balance Trade-off")
  legend("topright", legend = unique(metrics$method),
         col = 1:length(unique(metrics$method)), pch = 19)

  dev.off()

  # Plot 5: Coverage by sample size for each CI method
  pdf(file.path(output_dir, "coverage_by_n_and_ci.pdf"), width = 12, height = 10)

  par(mfrow = c(2, 2))
  for (i in seq_along(ci_methods)) {
    cov_col <- paste0("coverage_", ci_methods[i])
    if (cov_col %in% names(metrics)) {
      # Aggregate by n and method
      cov_by_n <- aggregate(
        as.formula(paste(cov_col, "~ n + method")),
        data = metrics,
        FUN = mean
      )

      # Create plot
      plot(1, type = "n", xlim = range(cov_by_n$n), ylim = c(0.8, 1),
           xlab = "Sample Size", ylab = "Coverage",
           main = paste("Coverage by Sample Size -", ci_labels[i]))

      for (m in unique(cov_by_n$method)) {
        subset_data <- cov_by_n[cov_by_n$method == m, ]
        lines(subset_data$n, subset_data[[cov_col]],
              col = which(unique(cov_by_n$method) == m),
              type = "b", pch = 19)
      }

      abline(h = 0.95, lty = 2, col = "red")
      if (i == 1) {
        legend("bottomright", legend = unique(cov_by_n$method),
               col = 1:length(unique(cov_by_n$method)),
               lty = 1, pch = 19, cex = 0.8)
      }
    }
  }

  dev.off()

  # Plot 6: Bootstrap convergence diagnostics
  pdf(file.path(output_dir, "bootstrap_convergence.pdf"), width = 10, height = 6)

  boot_success <- aggregate(
    n_boot_valid ~ method + n,
    data = results,
    FUN = function(x) mean(x == 500, na.rm = TRUE)
  )

  barplot(boot_success$n_boot_valid,
          names.arg = paste(boot_success$method, boot_success$n, sep = "-"),
          main = "Proportion of Successful Bootstrap Iterations",
          xlab = "Method-Sample Size", ylab = "Proportion",
          col = rep(rainbow(8), each = 4),
          las = 2)
  abline(h = 1, lty = 2, col = "green")

  dev.off()
}

# =============================================================================
# SECTION 9: MAIN EXECUTION SCRIPT
# =============================================================================

#' Main function to run complete simulation study
#'
#' @description Master function that orchestrates the entire simulation study
#' from parameter setup through results generation and reporting.
#'
#' @param config_file Character. Optional path to configuration file that
#'   sources custom simulation parameters
#'
#' @return List containing:
#'   \item{results}{Raw simulation results data frame}
#'   \item{metrics}{Aggregated performance metrics}
#'
#' @details Workflow:
#' 1. Load configuration (if provided)
#' 2. Run full simulation study with parallel processing
#' 3. Calculate performance metrics
#' 4. Generate summary tables
#' 5. Create diagnostic plots
#' 6. Save session info for reproducibility
#'
#' Runtime depends on:
#' - Number of scenarios (factorial design)
#' - Number of replications (n_sim)
#' - Number of bootstrap iterations (n_boot)
#' - Number of CPU cores (n_cores)
#'
#' @examples
#' # Run with default parameters
#' results <- main()
#'
#' # Run with custom configuration
#' results <- main("custom_sim_params.R")
#'
#' # Access results
#' head(results$metrics)
#' summary(results$results$hr_est)
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
#' @description Formats key simulation results as LaTeX tables suitable for
#' direct inclusion in manuscripts.
#'
#' @param metrics Data frame. Performance metrics from simulation
#' @param output_file Character. LaTeX file name (default "simulation_results.tex")
#'
#' @details Creates formatted tables for:
#' - Bias by method and sample size
#' - Coverage probability by method and CI type
#' - Other key performance metrics
#'
#' Tables use booktabs package for professional formatting.
#'
#' @examples
#' metrics <- calculate_performance_metrics(results)
#' export_latex_tables(metrics, "manuscript_tables.tex")
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
#' @description Performs quality checks on simulation results to identify
#' potential issues or unexpected patterns.
#'
#' @param results Data frame. Raw simulation results
#'
#' @return List of validation checks with pass/fail status
#'
#' @details Checks include:
#' - Completeness: All methods present for each scenario
#' - Missing values: No missing HR estimates
#' - Range checks: HR estimates in reasonable range
#' - Constraint checks: ESS ≤ n, positive weights
#'
#' @examples
#' validation <- validate_results(simulation_results)
#'
#' # Check if all validations passed
#' all(unlist(validation))
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
