# =============================================================================
# MAIC Simulation Study for Time-to-Event Outcomes with IPTW ATUT Comparison
#
# This code implements the simulation study described in the manuscript:
# "Comparative Evaluation of Weight Normalization Methods in Matching-Adjusted
# Indirect Comparison for Time-to-Event Outcomes"
# WITH IPTW ATUT (Average Treatment Effect in the Untreated) COMPARISON
#
# Author: Gregory Chen, etc
# Date: Aug 17, 2025
# R Version: 4.5.1
# =============================================================================
# KEY CHANGES IN THIS VERSION
# =============================================================================
#
# 1. IPTW NOW TARGETS ATUT (not ATT):
#    - ATUT = Average Treatment Effect in the Untreated (external/control population)
#    - This matches MAIC's target estimand
#    - External group (control, coded 0): weight = 1
#    - Internal group (treatment, coded 1): weight = (1-e)/e where e = P(Internal|X)
#
# 2. TRUE MARGINAL ATUT HR:
#    - Replaced conditional HR with marginal ATUT HR as the truth
#    - Computed via Monte Carlo simulation accounting for:
#      * Non-collapsibility of HR
#      * Censoring and dropout patterns
#      * Maximum follow-up time
#    - Used for bias and coverage evaluation
#
# =============================================================================
# MATHEMATICAL NOTATION ALIGNMENT WITH MANUSCRIPT
# =============================================================================
#
# Manuscript Notation -> Code Variable Mapping:
# ---------------------------------------------
# μ*_X2 (external means)     -> X_external_means, colMeans(X_external)
# X_i (covariates)          -> X_internal, X_external matrices
# ω(X_i) (weights)          -> weights, maic_weights, iptw_weights
# α (weight parameters)      -> alpha in optimization functions
# β (prognostic effects)     -> beta vector (β1, β2, β3, 0)
# γ (treatment effect)       -> gamma = log(conditional_hr)
# Δ (ATUT estimand)         -> true_marginal_hr (target of estimation)
# n (sample size)           -> n per arm
# ESS (effective sample)     -> ess = (Σw_i)²/Σw_i²
# d (covariate shift)       -> d parameter controlling population difference
#
# =============================================================================
# HIGH-LEVEL DESIGN SUMMARY
# =============================================================================
#
# TARGET ESTIMAND: ATUT (Average Treatment Effect in the Untreated)
# - Target population: EXTERNAL group (control)
# - Both MAIC and IPTW weight INTERNAL group to match EXTERNAL population
# - INTERNAL group = Treatment/Investigational (coded as 1)
# - EXTERNAL group = Control (coded as 0)
#
# WEIGHTING SCHEMES:
# - MAIC: Weights for INTERNAL group via method of moments to match EXTERNAL means
# - MEW: Maximum ESS weights for INTERNAL group
# - IPTW ATUT: Weights = (1-e)/e for INTERNAL, where e = P(Internal|X)
# - All methods: EXTERNAL group receives weight = 1
#
# TRUE PARAMETER:
# - Marginal ATUT HR computed via Monte Carlo
# - Accounts for non-collapsibility and censoring
# - Used as truth for bias and coverage evaluation
#
# =============================================================================
# ARCHITECTURE AND FUNCTION HIERARCHY
# =============================================================================
#
# 1. DATA GENERATION LAYER
# -------------------------
# Core functions for generating simulation data:
#
# generate_covariates(n, population, d)
#   ├── Used by: run_single_simulation_with_params(), calibrate_d_for_ess()
#   └── Purpose: Creates baseline covariates X1-X4 with controlled imbalance
#
# generate_survival_times(n, X, beta, gamma, Z, distribution, base_params)
#   ├── Used by: run_single_simulation_with_params(), compute_true_marginal_atut_hr()
#   └── Purpose: Generates event times from Weibull or log-normal distributions
#
# generate_censoring(n, T, dropout_rate, max_followup)
#   ├── Used by: run_single_simulation_with_params(), compute_true_marginal_atut_hr()
#   └── Purpose: Creates realistic dropout patterns (60% early, 40% late)
#
# 2. CALIBRATION LAYER
# --------------------
# Functions for parameter calibration:
#
# calibrate_d_for_ess(target_ess_ratio, n_calib, tol, method)
#   ├── Used by: run_simulation_study()
#   ├── Purpose: Finds d value to achieve target ESS ratio
#   └── Methods: bisection or Newton-Raphson
#
# calibrate_base_params(distribution, median_surv, surv_target, target_time)
#   ├── Used by: run_single_simulation_with_params(), compute_true_marginal_atut_hr()
#   └── Purpose: Sets Weibull/log-normal parameters for target survival
#
# calibrate_beta(target_event_rate, n_calib, distribution, base_params, ...)
#   ├── Used by: run_simulation_study()
#   └── Purpose: Calibrates prognostic effects for target event rate
#
# compute_true_marginal_atut_hr(conditional_hr, d, beta, distribution, ...)
#   ├── Used by: run_simulation_study()
#   ├── Purpose: Computes true marginal HR via Monte Carlo (10,000 samples)
#   └── Key feature: Accounts for non-collapsibility
#
# 3. WEIGHT CALCULATION LAYER
# ---------------------------
# Functions for computing different weight types:
#
# calculate_iptw_atut_weights(X_internal, X_external, X_for_ps, max_weight)
#   ├── Used by: run_single_simulation_with_params(), calibrate_d_for_ess()
#   ├── Purpose: Computes IPTW weights for ATUT: (1-e)/e
#   └── Note: e = P(Internal|X), weights internal to match external
#
# calculate_maic_weights(X_internal, X_external_means, max_weight)
#   ├── Used by: run_single_simulation_with_params()
#   ├── Purpose: Method of moments optimization for MAIC weights
#   └── Uses: BFGS optimization with gradient
#
# calculate_mew_weights(X_internal, X_external_means, max_weight)
#   ├── Used by: run_single_simulation_with_params()
#   ├── Purpose: Maximum ESS weights via quadratic programming
#   └── Uses: Closed-form solution with regularization
#
# normalize_weights(weights, method)
#   ├── Used by: run_single_simulation_with_params()
#   ├── Purpose: Applies normalization schemes
#   └── Methods: OW (original), SW1 (sum-to-1), SWN (sum-to-N), SWESS (sum-to-ESS)
#
# 4. STATISTICAL ANALYSIS LAYER
# -----------------------------
# Functions for outcome analysis:
#
# fit_weighted_cox(time, event, treatment, weights, n_boot, ...)
#   ├── Used by: run_single_simulation_with_params()
#   ├── Purpose: Cox regression with proper bootstrap
#   ├── Features: Recalculates weights in each bootstrap iteration
#   └── Returns: HR estimates with 4 CI types (robust, percentile, normal, BCa)
#
# calculate_ess(weights)
#   ├── Used by: run_single_simulation_with_params()
#   └── Purpose: Computes effective sample size
#
# calculate_balance(X_internal, X_external_means, weights)
#   ├── Used by: run_single_simulation_with_params()
#   └── Purpose: Assesses covariate balance via standardized differences
#
# 5. SIMULATION ORCHESTRATION LAYER
# ---------------------------------
# Main simulation control functions:
#
# run_single_simulation_with_params(params, sim_params)
#   ├── Used by: run_simulation_study() via parallel processing
#   ├── Purpose: Single Monte Carlo replication
#   ├── Process: Generate data → Calculate weights → Fit models → Return results
#   └── Returns: DataFrame with 12 methods × performance metrics
#
# run_simulation_study(sim_params, save_intermediate, save_dir)
#   ├── Used by: main()
#   ├── Purpose: Full factorial simulation with parallel processing
#   ├── Process:
#   │   1. Calibrate d values for ESS ratios
#   │   2. Create factorial design
#   │   3. Calibrate beta coefficients
#   │   4. Compute true marginal HRs
#   │   5. Run parallel simulations
#   │   6. Save intermediate results
#   └── Returns: Complete simulation results DataFrame
#
# 6. PERFORMANCE EVALUATION LAYER
# -------------------------------
# Functions for analyzing results:
#
# calculate_performance_metrics(results)
#   ├── Used by: main()
#   ├── Purpose: Computes bias, coverage, MSE using true marginal HR
#   └── Returns: Summary statistics by scenario and method
#
# create_comparison_tables(results)
#   ├── Used by: main()
#   ├── Purpose: MAIC vs IPTW ATUT comparisons
#   └── Returns: Overall comparison, best methods, relative efficiency
#
# create_summary_tables(results, output_format)
#   ├── Used by: main() or standalone
#   ├── Purpose: Formatted tables for reporting
#   └── Formats: LaTeX or standard DataFrame
#
# 7. DIAGNOSTIC AND HELPER FUNCTIONS
# ----------------------------------
# Analysis and debugging utilities:
#
# check_weight_distributions(results)
#   ├── Standalone use: check_weight_distributions(simulation_results)
#   ├── Purpose: Analyzes ESS and max weight by method
#   └── Returns: Weight statistics summary
#
# validate_bootstrap(results, expected_n_boot)
#   ├── Standalone use: validate_bootstrap(simulation_results, 500)
#   ├── Purpose: Checks bootstrap convergence
#   └── Returns: Bootstrap validity rates by scenario
#
# compare_methods(results)
#   ├── Standalone use: compare_methods(simulation_results)
#   ├── Purpose: Direct MAIC vs IPTW ATUT comparison
#   └── Returns: Comparative performance metrics
#
# summarize_convergence(results, n_sim)
#   ├── Standalone use: summarize_convergence(results, 10000)
#   ├── Purpose: Assesses Monte Carlo error
#   └── Returns: Convergence diagnostics with flags
#
# 8. MAIN ORCHESTRATOR
# -------------------
#
# main(sim_params, run_analysis)
#   ├── Entry point for full simulation study
#   ├── Process:
#   │   1. Validates parameters
#   │   2. Runs simulation via run_simulation_study()
#   │   3. Performs all analyses if run_analysis=TRUE
#   │   4. Saves results and summaries
#   └── Returns: List with results, analysis summaries, runtime
#
# =============================================================================

# Load required libraries
library(survival)
library(parallel)
library(MASS)  # For ginv if needed

# Set global parameters with flexible file naming
set.seed(12345)  # For reproducibility

# Enhanced simulation parameters with file naming options
sim_params <- list(
  # Sample sizes per arm
  n_sizes = c(20, 50, 200, 500),

  # Conditional HR values (exp(gamma) in the manuscript)
  conditional_hr = c(0.70, 0.85),  # Renamed from true_hr for clarity

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
  max_weight = 100,  # Maximum weight cap for stability
  n_mc_truth = 10000,  # Sample size for computing true marginal HR

  # Parallel computing
  n_cores = 4,  # Adjust based on your system

  # File naming options
  results_prefix = "maic_iptw_atut_sim",  # Updated prefix
  timestamp_files = TRUE  # Add timestamp to filenames
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Generate baseline covariates for internal and external populations
#'
#' @description Creates four baseline covariates (X1-X4) with controlled
#' differences between internal and external populations.
#'
#' Manuscript Reference: Section 3.1.1 "Baseline Covariate Generation"
#' - X1: Continuous covariate ~ N(0.5, 0.5²) for external
#' - X2-X4: Binary covariates ~ Bernoulli(0.5) for external
#' - Internal population: means scaled by parameter d
#' - d controls the degree of population difference (d=1 means identical)
#'
#' @param n Sample size to generate
#' @param population Either "external" (control) or "internal" (treatment)
#' @param d Scaling factor for population differences (d > 1 increases difference)
#' @return Matrix with columns X1, X2, X3, X4
generate_covariates <- function(n, population = "external", d = 1) {
  if (population == "external") {
    # External (control) population: Reference population
    # Manuscript: "X1 ~ N(0.5, 0.5²), X2-X4 ~ Bernoulli(0.5)"
    X1 <- rnorm(n, mean = 0.5, sd = 0.5)        # Continuous covariate
    X2 <- rbinom(n, size = 1, prob = 0.5)       # Binary prognostic factor
    X3 <- rbinom(n, size = 1, prob = 0.5)       # Binary prognostic factor
    X4 <- rbinom(n, size = 1, prob = 0.5)       # Binary non-prognostic factor
  } else {
    # Internal (treatment) population: means scaled by d
    # Manuscript: "introduce scaling factor d to control population differences"
    # As d increases from 1, populations become more different
    X1 <- rnorm(n, mean = 0.5 * d, sd = 0.5)    # Shifted mean
    X2 <- rbinom(n, size = 1, prob = min(0.5 * d, 1))  # Shifted probability
    X3 <- rbinom(n, size = 1, prob = min(0.5 * d, 1))  # Shifted probability
    X4 <- rbinom(n, size = 1, prob = min(0.5 * d, 1))  # Shifted probability
  }

  # Return as matrix for matrix operations in weight calculations
  return(cbind(X1 = X1, X2 = X2, X3 = X3, X4 = X4))
}

#' Calculate IPTW ATUT weights
#'
#' @description Computes inverse probability of treatment weights for ATUT
#'
#' Manuscript Reference: Section 2.7 "Relationship with IPTW"
#' For ATUT (Average Treatment Effect in the Untreated):
#' - Target population: EXTERNAL (control) group
#' - Weights for internal: ω(X) = (1-e(X))/e(X)
#' - Where e(X) = P(Internal|X) = P(Treatment=1|X)
#' - External group receives weight = 1 (reference population)
#'
#' Mathematical derivation:
#' ATUT requires E[Y(1)|External] which we estimate by weighting Internal to External
#' ω = P(External|X)/P(Internal|X) = (1-e(X))/e(X)
#'
#' @param X_internal Matrix of internal study covariates (treatment group)
#' @param X_external Matrix of external study covariates (control group)
#' @param X_for_ps Columns of X to use for propensity score model (default 1:3)
#' @param max_weight Maximum weight cap for numerical stability (default 100)
#' @return Vector of weights for internal population
calculate_iptw_atut_weights <- function(X_internal, X_external, X_for_ps = 1:3, max_weight = 100) {
  n_int <- nrow(X_internal)
  n_ext <- nrow(X_external)

  # Input validation
  if (n_int == 0 || n_ext == 0) {
    warning("Empty input to IPTW weight calculation")
    return(rep(1, n_int))  # Return unit weights as fallback
  }

  # Combine data for propensity score model
  # Manuscript: "logit P(T=2|X) = α₀ + α'X" where T=2 is external
  # We code: Internal = 1 (treatment), External = 0 (control)
  X_combined <- rbind(X_internal[, X_for_ps, drop = FALSE],
                      X_external[, X_for_ps, drop = FALSE])

  # Treatment indicator: 1 = internal (treatment), 0 = external (control)
  T_combined <- c(rep(1, n_int), rep(0, n_ext))

  # Create data frame for GLM
  # Handle both single and multiple covariates
  if (length(X_for_ps) == 1) {
    ps_data <- data.frame(
      T = T_combined,
      X1 = X_combined
    )
  } else {
    ps_data <- data.frame(
      T = T_combined,
      X_combined
    )
    colnames(ps_data) <- c("T", paste0("X", 1:length(X_for_ps)))
  }

  # Fit propensity score model
  # Manuscript: "propensity score e(X) = P(Internal|X)"
  ps_model <- tryCatch({
    glm(T ~ ., data = ps_data, family = binomial())
  }, warning = function(w) {
    # Suppress convergence warnings but still fit
    suppressWarnings(glm(T ~ ., data = ps_data, family = binomial()))
  }, error = function(e) {
    warning(paste("Propensity score model failed:", e$message))
    return(NULL)
  })

  if (is.null(ps_model)) {
    warning("Propensity score model failed, returning equal weights")
    return(rep(1, n_int))
  }

  # Get propensity scores for INTERNAL population
  # These are e(X) = P(Internal|X) for each internal observation
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
    return(rep(0.5, n_int))  # Default to 0.5 if prediction fails
  })

  # Calculate IPTW ATUT weights: (1-e)/e for internal
  # Manuscript equation: ω(X) = P(External|X)/P(Internal|X) = (1-e(X))/e(X)
  # This weights internal group to match external population distribution
  weights <- (1 - ps_internal) / ps_internal

  # Handle extreme weights for numerical stability
  # Manuscript: "cap extreme weights at max_weight"
  if (any(is.na(weights)) || any(is.infinite(weights))) {
    warning("NA or infinite weights detected, setting to 1")
    weights[is.na(weights) | is.infinite(weights)] <- 1
  }

  # Weight truncation to prevent single observations dominating
  weights <- pmin(weights, max_weight)

  return(weights)
}

#' Compute true marginal ATUT HR via Monte Carlo
#'
#' @description Computes the true marginal HR for ATUT estimand
#' accounting for non-collapsibility and censoring patterns
#'
#' Manuscript Reference: Section 2.1 "Target Estimand"
#' The marginal ATUT HR differs from conditional HR due to non-collapsibility
#' of the hazard ratio when there are prognostic covariates.
#'
#' Algorithm:
#' 1. Generate large sample from EXTERNAL population
#' 2. Randomly assign to treatment/control (ensures same covariate distribution)
#' 3. Generate survival times with treatment effect
#' 4. Apply realistic censoring
#' 5. Fit unweighted Cox model to get marginal HR
#'
#' This is the TRUE parameter that methods are trying to estimate
#'
#' @param conditional_hr Conditional HR (exp(gamma) in Cox model)
#' @param d Covariate shift parameter between populations
#' @param beta Vector of prognostic effects (β₁, β₂, β₃, 0)
#' @param distribution "weibull" or "lognormal"
#' @param base_params Baseline hazard parameters (shape, scale or mu, sigma)
#' @param dropout_rate Proportion of random dropout
#' @param max_followup Administrative censoring time
#' @param n_mc Monte Carlo sample size (default 10000)
#' @return True marginal ATUT HR
compute_true_marginal_atut_hr <- function(conditional_hr, d, beta, distribution,
                                          base_params, dropout_rate, max_followup = 48,
                                          n_mc = 10000) {
  # Set seed for reproducibility within this function
  # Ensures same true HR for same parameter combination
  set.seed(123456)

  # Generate large sample from EXTERNAL population
  # This is our target population for ATUT
  X_external_mc <- generate_covariates(n_mc, "external")

  # Randomly assign half to each treatment group
  # Key insight: Both groups now have EXTERNAL covariate distribution
  # This mimics what perfect weighting would achieve
  n_half <- n_mc / 2
  idx_treat <- sample(1:n_mc, n_half)
  idx_control <- setdiff(1:n_mc, idx_treat)

  X_treat <- X_external_mc[idx_treat, ]
  X_control <- X_external_mc[idx_control, ]

  # Calculate gamma from conditional HR
  # Manuscript: γ = log(HR_conditional)
  gamma <- log(conditional_hr)

  # Generate survival times for treatment group (with treatment effect)
  # Manuscript equation: h(t|X,Z) = h₀(t)exp(β'X + γZ)
  T_treat <- generate_survival_times(n_half, X_treat, beta, gamma, rep(1, n_half),
                                     distribution, base_params)

  # Generate survival times for control group (no treatment effect)
  T_control <- generate_survival_times(n_half, X_control, beta, 0, rep(0, n_half),
                                       distribution, base_params)

  # Apply censoring to mimic real trial conditions
  # Important: censoring pattern affects marginal HR
  C_treat <- generate_censoring(n_half, T_treat, dropout_rate, max_followup)
  C_control <- generate_censoring(n_half, T_control, dropout_rate, max_followup)

  # Observed outcomes
  Y_treat <- pmin(T_treat, C_treat)
  Y_control <- pmin(T_control, C_control)
  delta_treat <- as.numeric(T_treat <= C_treat)
  delta_control <- as.numeric(T_control <= C_control)

  # Combine data
  Y_combined <- c(Y_treat, Y_control)
  delta_combined <- c(delta_treat, delta_control)
  Z_combined <- c(rep(1, n_half), rep(0, n_half))

  # Fit unweighted Cox model to get marginal HR
  # This is the parameter our methods are trying to estimate
  surv_obj <- Surv(Y_combined, delta_combined)
  cox_fit <- tryCatch({
    coxph(surv_obj ~ Z_combined)
  }, error = function(e) {
    warning(paste("Cox model failed in true HR computation:", e$message))
    return(NULL)
  })

  if (is.null(cox_fit)) {
    warning("Using conditional HR as fallback for true marginal HR")
    return(conditional_hr)
  }

  # Extract marginal HR
  # This accounts for non-collapsibility: marginal HR ≠ conditional HR
  marginal_hr <- exp(coef(cox_fit))

  return(as.numeric(marginal_hr))
}

#' Calibrate d value to achieve target ESS ratio
calibrate_d_for_ess <- function(target_ess_ratio, n_calib = 10000, tol = 0.01,
                                method = c("bisection", "newton-ralphson")) {
  # Function to compute ESS ratio for given d
  compute_ess_ratio <- function(d) {
    # Generate large samples for calibration
    X_external <- generate_covariates(n_calib, "external")
    X_internal <- generate_covariates(n_calib, "internal", d)

    # Compute IPTW ATUT weights
    weights <- calculate_iptw_atut_weights(X_internal, X_external, X_for_ps = 1:3)

    # Weights for combined population (internal gets IPTW weights, external gets 1)
    weights_combined <- c(weights, rep(1, n_calib))

    # Compute ESS
    ess <- sum(weights_combined)^2 / sum(weights_combined^2)

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

      if (ess_mid > target_ess_ratio) {
        d_lower <- d_mid
      } else {
        d_upper <- d_mid
      }
    }

    d_current <- (d_lower + d_upper) / 2
  } else if(method == "newton-ralphson") {
    # Newton-Raphson method
    d_current <- 1.1
    max_iter <- 50
    h <- 0.001

    ess_at_1 <- compute_ess_ratio(1.0)
    if (ess_at_1 < target_ess_ratio) {
      warning(sprintf("Target ESS ratio %.3f is higher than maximum achievable %.3f at d=1",
                      target_ess_ratio, ess_at_1))
      return(1.0)
    }

    for (iter in 1:max_iter) {
      f_current <- compute_ess_ratio(d_current) - target_ess_ratio

      if (abs(f_current) < tol) {
        break
      }

      f_plus <- compute_ess_ratio(d_current + h) - target_ess_ratio
      f_minus <- compute_ess_ratio(d_current - h) - target_ess_ratio
      f_prime <- (f_plus - f_minus) / (2 * h)

      if (abs(f_prime) < 1e-10) {
        if (f_current > 0) {
          d_current <- d_current + 0.1
        } else {
          d_current <- d_current - 0.05
        }
      } else {
        d_new <- d_current - f_current / f_prime
        d_new <- pmax(1.001, pmin(2.0, d_new))

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
    # Ensure valid parameters
    if (base_params$shape <= 0 || base_params$lambda <= 0) {
      stop(sprintf("Invalid Weibull parameters: shape=%.4f, lambda=%.4f",
                   base_params$shape, base_params$lambda))
    }

    # Generate from Weibull with covariates
    U <- runif(n)
    # T = (-log(U) / (lambda * exp(lp)))^(1/shape)
    T <- (-log(U) / (base_params$lambda * exp(lp)))^(1/base_params$shape)

  } else if (distribution == "lognormal") {
    # Generate from log-normal with covariates
    # log(T) ~ N(mu - lp, sigma^2)
    T <- exp(rnorm(n, mean = base_params$mu - lp, sd = base_params$sigma))
  } else {
    stop("Unknown distribution")
  }

  # Ensure positive survival times
  T <- pmax(T, 0.001)

  return(as.numeric(T))
}

#' Calibrate baseline distribution parameters
calibrate_base_params <- function(distribution, median_surv, surv_target, target_time = 48) {
  if (distribution == "weibull") {
    # Ensure valid survival probability
    if (surv_target <= 0 || surv_target >= 1) {
      surv_target <- pmax(0.001, pmin(0.999, surv_target))
    }

    # For Weibull: S(t) = exp(-(lambda * t)^shape)
    # At median: S(median_surv) = 0.5
    # At target_time: S(target_time) = surv_target

    # From S(median_surv) = 0.5: -log(0.5) = (lambda * median_surv)^shape
    # From S(target_time) = surv_target: -log(surv_target) = (lambda * target_time)^shape

    # Taking ratio and solving for shape:
    # log(-log(0.5)) / log(-log(surv_target)) = log(median_surv/target_time) * shape / log(target_time/target_time)

    # Corrected calculation:
    shape <- log(-log(surv_target)) / log(target_time / median_surv)

    # Ensure shape is positive
    if (shape <= 0) {
      warning(sprintf("Invalid Weibull shape (%.4f). Using shape = 1 (exponential).", shape))
      shape <- 1
    }

    # Calculate lambda from median survival
    lambda <- (-log(0.5))^(1/shape) / median_surv

    return(list(lambda = lambda, shape = shape))

  } else if (distribution == "lognormal") {
    # For log-normal: log(T) ~ N(mu, sigma^2)
    # Median survival: exp(mu) = median_surv, so mu = log(median_surv)
    mu <- log(median_surv)

    # For survival probability at target_time:
    # P(T > target_time) = P(log(T) > log(target_time)) = 1 - Phi((log(target_time) - mu) / sigma)
    # surv_target = 1 - Phi((log(target_time) - mu) / sigma)
    # Phi^(-1)(1 - surv_target) = (log(target_time) - mu) / sigma
    sigma <- (log(target_time) - mu) / qnorm(1 - surv_target)

    # Ensure sigma is positive
    if (sigma <= 0) {
      warning(sprintf("Invalid log-normal sigma (%.4f). Using sigma = 0.5.", sigma))
      sigma <- 0.5
    }

    return(list(mu = mu, sigma = sigma))
  }
}



#' Generate censoring times based on dropout pattern
generate_censoring <- function(n, T, dropout_rate, max_followup = 48) {
  n_dropout <- rbinom(1, n, dropout_rate)

  if (n_dropout == 0) {
    return(rep(max_followup, n))
  }

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

  return(C)
}

# Also fix the calibrate_beta function to handle extreme cases better
calibrate_beta <- function(target_event_rate, n_calib = 5000, distribution,
                           base_params, dropout_rate, max_followup = 48,
                           tol = 0.01, method = c("bisection", "newton-raphson")) {

  compute_event_rate <- function(beta_scale) {
    beta <- c(beta_scale, beta_scale, beta_scale, 0)

    X <- generate_covariates(n_calib, "external")
    T <- generate_survival_times(n_calib, X, beta, gamma = 0, Z = rep(0, n_calib),
                                 distribution, base_params)

    C <- generate_censoring(n_calib, T, dropout_rate, max_followup)

    Y <- pmin(T, C)
    delta <- as.numeric(T <= C)

    event_rate <- mean(delta)
    return(event_rate)
  }

  method <- match.arg(method)

  # Start with reasonable bounds
  beta_lower <- -3
  beta_upper <- 3

  # Test the bounds
  rate_lower <- compute_event_rate(beta_lower)
  rate_upper <- compute_event_rate(beta_upper)

  # Adjust bounds if necessary
  max_attempts <- 5
  attempt <- 0

  while (attempt < max_attempts && (rate_lower > target_event_rate || rate_upper < target_event_rate)) {
    attempt <- attempt + 1

    if (rate_lower > target_event_rate) {
      # Need more negative beta
      beta_lower <- beta_lower - 1
      rate_lower <- compute_event_rate(beta_lower)
    }

    if (rate_upper < target_event_rate) {
      # Need more positive beta
      beta_upper <- beta_upper + 1
      rate_upper <- compute_event_rate(beta_upper)
    }
  }

  # If still can't bracket, use closest bound
  if (target_event_rate < rate_lower) {
    return(c(beta_lower, beta_lower, beta_lower, 0))
  }
  if (target_event_rate > rate_upper) {
    return(c(beta_upper, beta_upper, beta_upper, 0))
  }

  if(method == "bisection") {
    # Bisection method
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

  } else if(method == "newton-raphson") {
    # Newton-Raphson method
    beta_current <- 0.0
    max_iter <- 50
    h <- 0.01

    for (iter in 1:max_iter) {
      f_current <- compute_event_rate(beta_current) - target_event_rate

      if (abs(f_current) < tol) {
        break
      }

      # Numerical derivative
      rate_plus <- compute_event_rate(beta_current + h)
      rate_minus <- compute_event_rate(beta_current - h)
      f_prime <- (rate_plus - rate_minus) / (2 * h)

      if (abs(f_prime) < 1e-10) {
        # Derivative too small, use bisection fallback
        if (f_current > 0) {
          beta_current <- beta_current - 0.1
        } else {
          beta_current <- beta_current + 0.1
        }
      } else {
        beta_new <- beta_current - f_current / f_prime
        # Constrain to reasonable range
        beta_new <- pmax(beta_lower, pmin(beta_upper, beta_new))
        beta_current <- beta_new
      }
    }

    beta_scale <- beta_current
  }

  return(c(beta_scale, beta_scale, beta_scale, 0))
}

#' Calculate MAIC weights using method of moments
#'
#' @description Computes weights via method of moments optimization
#' to balance covariate means between populations
#'
#' Manuscript Reference: Section 2.3 "MAIC Weights"
#' Solves: min Σ exp(α'X_i^c) subject to Σ w_i X_i = n × μ_external
#' where X_i^c = X_i - μ_external (centered covariates)
#'
#' The optimization problem is strictly convex with unique solution.
#' Weights have form: w_i = exp(α'X_i^c)
#'
#' @param X_internal Matrix of internal study covariates
#' @param X_external_means Vector of external population means
#' @param max_weight Maximum weight truncation (default 100)
#' @return Vector of weights for internal population
calculate_maic_weights <- function(X_internal, X_external_means, max_weight = 100) {
  n <- nrow(X_internal)
  p <- ncol(X_internal)

  # Center covariates around external means
  # Manuscript: X_i^c = X_i - μ*_X2
  X_centered <- sweep(X_internal, 2, X_external_means, "-")

  # Method of moments optimization
  # Manuscript: Q(α) = Σ exp(α'X_i^c)
  objective <- function(alpha) {
    # Calculate log weights with numerical stability
    log_weights <- as.numeric(X_centered %*% alpha)
    log_weights <- pmin(log_weights, 20)  # Prevent exp overflow
    weights <- exp(log_weights)

    # Handle degenerate case
    if (sum(weights) == 0) return(1e10)

    # Check balance constraint: weighted means should equal external means
    weighted_means <- colSums(weights * X_internal) / sum(weights)

    # Return squared deviation from target (for minimization)
    sum((weighted_means - X_external_means)^2)
  }

  # Gradient of objective function
  # Manuscript: ∇Q(α) = Σ X_i^c exp(α'X_i^c)
  gradient <- function(alpha) {
    log_weights <- as.numeric(X_centered %*% alpha)
    log_weights <- pmin(log_weights, 20)
    weights <- exp(log_weights)

    if (sum(weights) == 0) return(rep(0, p))

    W <- sum(weights)
    weighted_X <- colSums(weights * X_internal)

    # Compute gradient components
    grad <- numeric(p)
    for (j in 1:p) {
      # Derivative of balance constraint violation
      term1 <- sum(weights * X_centered[, j] * X_internal[, j]) / W
      term2 <- weighted_X[j] * sum(weights * X_centered[, j]) / W^2
      grad[j] <- 2 * (weighted_X[j]/W - X_external_means[j]) * (term1 - term2)
    }
    return(grad)
  }

  # Optimize using BFGS (quasi-Newton method)
  # Starting from α = 0 (equal weights)
  result <- tryCatch({
    optim(par = rep(0, p),        # Initial α = 0
          fn = objective,          # Objective function
          gr = gradient,           # Gradient
          method = "BFGS",         # Quasi-Newton
          control = list(maxit = 1000))
  }, error = function(e) {
    # Fallback without gradient if analytical gradient fails
    optim(par = rep(0, p),
          fn = objective,
          method = "BFGS",
          control = list(maxit = 1000))
  })

  # Extract optimized weights
  log_weights <- as.numeric(X_centered %*% result$par)
  log_weights <- pmin(log_weights, 20)
  weights <- exp(log_weights)

  # Cap extreme weights for numerical stability
  # Manuscript: "truncate weights at max_weight"
  weights <- pmin(weights, max_weight)

  return(weights)
}

#' Calculate maximum ESS weights (MEW)
#'
#' @description Computes weights that maximize effective sample size
#' while achieving exact balance
#'
#' Manuscript Reference: Section 2.3.2 "Maximum ESS Weights"
#' Optimization: max ESS(w) = (Σw_i)²/Σw_i²
#' subject to: Σw_i = n, Σw_i X_i = n × μ_external
#'
#' Closed form solution: w = A(A'A)^(-1)b
#' where A = [1, X], b = [n, n×μ_external]
#'
#' @param X_internal Matrix of internal study covariates
#' @param X_external_means Vector of external population means
#' @param max_weight Maximum weight truncation (default 100)
#' @return Vector of weights for internal population
calculate_mew_weights <- function(X_internal, X_external_means, max_weight = 100) {
  n <- nrow(X_internal)
  p <- ncol(X_internal)

  # Construct constraint matrix
  # A = [1_n | X] for constraints on sum and means
  A <- cbind(rep(1, n), X_internal)

  # Right-hand side of constraints
  # b = [n, n×μ_external] for sum=n and means=μ_external
  b <- c(n, n * X_external_means)

  # Compute A'A for closed-form solution
  AtA <- t(A) %*% A

  # Check for near-singularity and add regularization if needed
  # This prevents numerical instability
  if (det(AtA) < 1e-10 || rcond(AtA) < 1e-10) {
    # Add small ridge penalty for numerical stability
    # Manuscript: "regularization parameter ε = 10^-6"
    AtA <- AtA + diag(1e-6, ncol(AtA))
  }

  # Closed-form solution for MEW
  # Manuscript: w = A(A'A)^(-1)b
  weights <- tryCatch({
    A %*% solve(AtA) %*% b
  }, error = function(e) {
    # Use generalized inverse if regular inverse fails
    A %*% MASS::ginv(AtA) %*% b
  })

  weights <- as.numeric(weights)

  # Ensure non-negative weights
  weights <- pmax(weights, 0)

  # Handle degenerate case
  if (sum(weights) == 0) {
    weights <- rep(1, n)
  }

  # Cap extreme weights
  weights <- pmin(weights, max_weight)

  return(weights)
}

#' Normalize weights according to different schemes
#'
#' @description Applies various normalization schemes to base weights
#'
#' Manuscript Reference: Section 2.4 "Weight Normalization Schemes"
#' Four schemes evaluated:
#' - OW: Original (unnormalized) weights
#' - SW1: Sum to one (Σw_i = 1)
#' - SWN: Sum to N (Σw_i = n)
#' - SWESS: Sum to ESS (Σw_i = ESS(w))
#'
#' @param weights Vector of base weights
#' @param method Normalization method: "OW", "SW1", "SWN", or "SWESS"
#' @return Normalized weight vector
normalize_weights <- function(weights, method) {
  n <- length(weights)

  # Handle edge case of zero weights
  if (sum(weights) == 0) {
    return(rep(1, n))
  }

  switch(method,
         # Original weights (no normalization)
         "OW" = weights,

         # Sum to one (standard normalization)
         # Manuscript: w̃_i = w_i / Σw_j
         "SW1" = weights / sum(weights),

         # Sum to N (preserve sample size interpretation)
         # Manuscript: w̃_i = n × w_i / Σw_j
         "SWN" = weights * n / sum(weights),

         # Sum to ESS (balance between SW1 and SWN)
         # Manuscript: w̃_i = ESS(w) × w_i / Σw_j
         "SWESS" = {
           ess <- sum(weights)^2 / sum(weights^2)
           weights * ess / sum(weights)
         },

         stop("Unknown normalization method")
  )
}

#' Calculate effective sample size
#'
#' @description Computes ESS as measure of weight variability
#'
#' Manuscript Reference: Section 2.5 "Effective Sample Size"
#' ESS = (Σw_i)²/Σw_i²
#'
#' Interpretation:
#' - ESS = n: uniform weights (no loss of precision)
#' - ESS < n: variable weights (loss of precision)
#' - ESS/n: efficiency relative to unweighted analysis
#'
#' @param weights Vector of weights
#' @return Effective sample size
calculate_ess <- function(weights) {
  if (sum(weights^2) == 0) return(0)

  # Standard ESS formula
  # Manuscript equation (14)
  sum(weights)^2 / sum(weights^2)
}

#' Fit weighted Cox model with proper bootstrap
#'
#' @description Fits Cox regression with various CI methods
#' Bootstrap recalculates weights in each iteration for valid inference
#'
#' Manuscript Reference: Section 2.6 "Variance Estimation"
#' Four CI methods implemented:
#' 1. Robust sandwich estimator
#' 2. Bootstrap percentile
#' 3. Bootstrap normal approximation
#' 4. BCa (bias-corrected and accelerated) bootstrap
#'
#' Key feature: Weights are recalculated in each bootstrap sample
#' to properly account for weight estimation uncertainty
#'
#' @param time Observed survival times
#' @param event Event indicators (1=event, 0=censored)
#' @param treatment Treatment indicators (1=internal, 0=external)
#' @param weights Analysis weights (NULL for unweighted)
#' @param n_boot Number of bootstrap iterations
#' @param X_internal Internal study covariates (for bootstrap)
#' @param X_external External study covariates (for bootstrap)
#' @param weight_type "maic", "mew", or "iptw"
#' @param weight_params List of weight calculation parameters
#' @return List with HR estimates and confidence intervals
fit_weighted_cox <- function(time, event, treatment, weights = NULL,
                             n_boot = 500, X_internal = NULL, X_external = NULL,
                             weight_type = "maic", weight_params = list()) {

  # Create survival object
  surv_obj <- Surv(time, event)

  # Initial Cox model fit
  cox_fit <- tryCatch({
    if (is.null(weights)) {
      coxph(surv_obj ~ treatment)
    } else {
      coxph(surv_obj ~ treatment, weights = weights, robust = TRUE)
    }
  }, error = function(e) {
    return(NULL)
  })

  # Handle fitting failure
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

  # Extract point estimate and robust SE
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

  # Robust confidence interval
  ci_robust_lower <- exp(log_hr - 1.96 * robust_se)
  ci_robust_upper <- exp(log_hr + 1.96 * robust_se)

  # Initialize bootstrap results
  ci_boot_perc_lower <- ci_boot_perc_upper <- NA
  ci_boot_norm_lower <- ci_boot_norm_upper <- NA
  ci_boot_bca_lower <- ci_boot_bca_upper <- NA
  boot_log_hr_dist <- NULL
  n_boot_valid <- 0

  # Bootstrap with weight recalculation
  if (n_boot > 0 && !is.null(X_internal) && !is.null(X_external)) {
    # Identify treatment groups
    trt_1_idx <- which(treatment == 1)  # Internal (treatment)
    trt_0_idx <- which(treatment == 0)  # External (control)
    n_1 <- length(trt_1_idx)
    n_0 <- length(trt_0_idx)

    # Check for valid groups
    if (n_0 > 0 && n_1 > 0) {
      boot_log_hr <- numeric(n_boot)
      boot_valid <- logical(n_boot)

      for (b in 1:n_boot) {
        # Stratified bootstrap sampling
        # Sample indices within each group
        boot_internal_idx <- sample(1:n_1, n_1, replace = TRUE)
        boot_external_idx <- sample(1:n_0, n_0, replace = TRUE)

        # Get the actual data indices
        boot_idx_1 <- trt_1_idx[boot_internal_idx]
        boot_idx_0 <- trt_0_idx[boot_external_idx]

        # Bootstrap covariates using the correct indices
        X_int_boot <- X_internal[boot_internal_idx, , drop = FALSE]
        X_ext_boot <- X_external[boot_external_idx, , drop = FALSE]

        # Recalculate weights for bootstrap sample
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
          rep(1, n_1)  # Fallback to unit weights
        })

        # Combine bootstrap data
        boot_time <- c(time[boot_idx_1], time[boot_idx_0])
        boot_event <- c(event[boot_idx_1], event[boot_idx_0])
        boot_treatment <- c(rep(1, n_1), rep(0, n_0))

        # Apply weights
        boot_weights_combined <- c(boot_weights, rep(1, n_0))

        # Fit Cox model on bootstrap sample
        boot_cox <- tryCatch({
          boot_surv <- Surv(boot_time, boot_event)
          coxph(boot_surv ~ boot_treatment,
                weights = boot_weights_combined,
                robust = TRUE)
        }, error = function(e) {
          return(NULL)
        })

        # Store bootstrap estimate if valid
        if (!is.null(boot_cox) && !is.na(coef(boot_cox))) {
          boot_log_hr[b] <- as.numeric(coef(boot_cox))
          boot_valid[b] <- TRUE
        } else {
          boot_valid[b] <- FALSE
        }
      }

      # Process bootstrap results
      boot_log_hr_valid <- boot_log_hr[boot_valid]
      n_boot_valid <- sum(boot_valid)

      # Calculate bootstrap CIs if sufficient valid samples
      if (n_boot_valid >= 20) {
        # Bootstrap percentile CI
        boot_hr <- exp(boot_log_hr_valid)
        ci_boot_perc_lower <- as.numeric(quantile(boot_hr, 0.025, na.rm = TRUE))
        ci_boot_perc_upper <- as.numeric(quantile(boot_hr, 0.975, na.rm = TRUE))

        # Bootstrap normal CI
        boot_se <- sd(boot_log_hr_valid, na.rm = TRUE)
        ci_boot_norm_lower <- exp(log_hr - 1.96 * boot_se)
        ci_boot_norm_upper <- exp(log_hr + 1.96 * boot_se)

        # BCa bootstrap CI (simplified)
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

  return(list(
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
  ))
}
#' Calculate balance metrics
#'
#' @description Assesses covariate balance after weighting
#'
#' Manuscript Reference: Section 2.8 "Balance Assessment"
#' Standardized difference: d_j = |μ_j,weighted - μ_j,external| / σ_j
#' where σ_j is SD in external population
#'
#' Threshold: |d_j| < 0.1 indicates adequate balance
#'
#' @param X_internal Internal study covariates
#' @param X_external_means External population means
#' @param weights Weights for internal population
#' @return Maximum standardized difference across covariates
calculate_balance <- function(X_internal, X_external_means, weights) {
  # Calculate weighted means for internal population
  weighted_means <- colSums(weights * X_internal) / sum(weights)

  # Standard deviations in external population (reference)
  # Manuscript: "σ for continuous, √(p(1-p)) for binary"
  # X1 is continuous with SD = 0.5
  # X2-X4 are binary with SD = √(0.5 × 0.5) = 0.5
  external_sds <- c(0.5, rep(sqrt(0.5 * (1 - 0.5)), ncol(X_internal) - 1))

  # Handle dimension mismatch
  if (length(external_sds) < length(weighted_means)) {
    external_sds <- rep(external_sds[1], length(weighted_means))
  }

  # Standardized differences
  # Manuscript equation: d_j = |μ_j,weighted - μ_j,external| / σ_j
  std_diffs <- abs(weighted_means - X_external_means) / external_sds

  # Return maximum imbalance
  return(max(std_diffs))
}

#' Run single simulation iteration
#'
#' @description Executes one Monte Carlo replication of the simulation
#'
#' Manuscript Reference: Section 3.2 "Simulation Procedure"
#' For each scenario:
#' 1. Generate data (covariates, survival times, censoring)
#' 2. Calculate weights (MAIC, MEW, IPTW)
#' 3. Apply normalizations (OW, SW1, SWN, SWESS)
#' 4. Fit Cox models with bootstrap
#' 5. Evaluate performance metrics
#'
#' Total of 12 methods tested:
#' - 4 MAIC normalizations
#' - 4 MEW normalizations
#' - 4 IPTW normalizations
#'
#' @param params List of scenario parameters
#' @param sim_params Global simulation parameters
#' @return DataFrame with results for all 12 methods
run_single_simulation_with_params <- function(params, sim_params) {
  # Extract scenario parameters
  n <- params$n                          # Sample size per arm
  conditional_hr <- params$conditional_hr # Conditional HR (exp(γ))
  true_marginal_hr <- params$true_marginal_hr # Pre-computed marginal ATUT HR
  distribution <- params$distribution     # Weibull or lognormal
  median_surv <- params$median_surv      # Median survival
  surv_48m <- params$surv_48m           # 48-month survival
  dropout_rate <- params$dropout_rate    # Dropout proportion
  event_rate <- params$event_rate        # Target event rate
  d <- params$d                          # Population difference parameter
  beta <- params$beta                    # Prognostic effects
  gamma <- params$gamma                  # log(conditional_hr)
  weighting_scenario <- params$weighting_scenario # correct/under/over
  max_followup <- params$max_followup    # Administrative censoring

  # Get simulation settings
  n_boot <- ifelse(is.null(sim_params$n_boot), 500, sim_params$n_boot)
  max_weight <- ifelse(is.null(sim_params$max_weight), 100, sim_params$max_weight)

  # ===========================================
  # STEP 1: Generate Data
  # ===========================================

  # Generate baseline distribution parameters
  base_params <- calibrate_base_params(distribution, median_surv, surv_48m,
                                       target_time = max_followup)

  # Generate covariates for both populations
  # Manuscript: "Internal = treatment, External = control"
  X_internal <- generate_covariates(n, "internal", d)
  X_external <- generate_covariates(n, "external")

  # Generate survival times
  # Internal: has treatment effect (γ)
  # External: no treatment effect (γ = 0)
  T_internal <- generate_survival_times(n, X_internal, beta, gamma, rep(1, n),
                                        distribution, base_params)
  T_external <- generate_survival_times(n, X_external, beta, 0, rep(0, n),
                                        distribution, base_params)

  # Apply censoring pattern
  C_internal <- generate_censoring(n, T_internal, dropout_rate, max_followup)
  C_external <- generate_censoring(n, T_external, dropout_rate, max_followup)

  # Create observed data
  Y_internal <- pmin(T_internal, C_internal)  # Observed time
  Y_external <- pmin(T_external, C_external)
  delta_internal <- as.numeric(T_internal <= C_internal)  # Event indicator
  delta_external <- as.numeric(T_external <= C_external)

  # Combine data for analysis
  # Treatment coding: Internal = 1, External = 0
  Y_combined <- c(Y_internal, Y_external)
  delta_combined <- c(delta_internal, delta_external)
  Z_combined <- c(rep(1, n), rep(0, n))

  # ===========================================
  # STEP 2: Select Covariates for Weighting
  # ===========================================

  # Determine which covariates to use based on scenario
  # Manuscript Section 3.1.4 "Model Specification Scenarios"
  if (weighting_scenario == "correct") {
    # Include all prognostic factors (X1, X2, X3)
    X_for_weights <- X_external[, 1:3]
    X_internal_for_weights <- X_internal[, 1:3]
    X_external_means <- colMeans(X_external[, 1:3])
    X_for_ps <- 1:3
  } else if (weighting_scenario == "under") {
    # Omit prognostic X3 (under-specified)
    X_for_weights <- X_external[, 1:2]
    X_internal_for_weights <- X_internal[, 1:2]
    X_external_means <- colMeans(X_external[, 1:2])
    X_for_ps <- 1:2
  } else if (weighting_scenario == "over") {
    # Include non-prognostic X4 (over-specified)
    X_for_weights <- X_external
    X_internal_for_weights <- X_internal
    X_external_means <- colMeans(X_external)
    X_for_ps <- 1:4
  }

  # ===========================================
  # STEP 3: Calculate Base Weights
  # ===========================================

  # MAIC weights (method of moments)
  maic_weights <- tryCatch({
    calculate_maic_weights(X_internal_for_weights, X_external_means, max_weight)
  }, error = function(e) {
    rep(1, n)  # Fallback to unit weights
  })

  # MEW weights (maximum ESS)
  mew_weights <- tryCatch({
    calculate_mew_weights(X_internal_for_weights, X_external_means, max_weight)
  }, error = function(e) {
    rep(1, n)
  })

  # IPTW ATUT weights
  iptw_weights <- tryCatch({
    calculate_iptw_atut_weights(X_internal, X_external, X_for_ps, max_weight)
  }, error = function(e) {
    rep(1, n)
  })

  # ===========================================
  # STEP 4: Initialize Results Storage
  # ===========================================

  results <- data.frame(
    n = integer(),
    conditional_hr = numeric(),
    true_marginal_hr = numeric(),  # True parameter for evaluation
    distribution = character(),
    median_surv = numeric(),
    surv_48m = numeric(),
    dropout_rate = numeric(),
    event_rate = numeric(),
    d = numeric(),
    weighting_scenario = character(),
    method = character(),
    hr_est = numeric(),             # Estimated HR
    ci_robust_lower = numeric(),    # Robust CI
    ci_robust_upper = numeric(),
    ci_boot_perc_lower = numeric(), # Bootstrap percentile CI
    ci_boot_perc_upper = numeric(),
    ci_boot_norm_lower = numeric(), # Bootstrap normal CI
    ci_boot_norm_upper = numeric(),
    ci_boot_bca_lower = numeric(),  # BCa CI
    ci_boot_bca_upper = numeric(),
    ess = numeric(),                # Effective sample size
    max_weight = numeric(),         # Maximum weight
    max_std_diff = numeric(),       # Balance metric
    n_boot_valid = numeric(),       # Valid bootstrap samples
    stringsAsFactors = FALSE
  )

  # ===========================================
  # STEP 5: Process MAIC and MEW Methods
  # ===========================================

  # Methods: OW, SW1, SWN, SWESS for both MAIC and MEW
  methods <- c("OW", "SW1", "SWN", "SWESS", "MEW", "MEW1", "MEWN", "MEWESS")

  for (method in methods) {
    # Determine base weights and normalization
    if (substr(method, 1, 3) == "MEW") {
      base_weights <- mew_weights
      weight_type <- "mew"
      # MEW -> OW, MEW1 -> SW1, MEWN -> SWN, MEWESS -> SWESS
      norm_method <- ifelse(method == "MEW", "OW",
                            ifelse(method == "MEW1", "SW1",
                                   ifelse(method == "MEWN", "SWN", "SWESS")))
    } else {
      base_weights <- maic_weights
      weight_type <- "maic"
      norm_method <- method
    }

    # Apply normalization to weights for INTERNAL group
    weights_internal <- normalize_weights(base_weights, norm_method)

    # Combined weights: Internal = weighted, External = 1 (reference)
    weights_combined <- c(weights_internal, rep(1, n))

    # Prepare parameters for Cox model fitting
    weight_params <- list(
      norm_method = norm_method,
      max_weight = max_weight,
      X_for_ps = X_for_ps
    )

    # Fit Cox model with bootstrap CI
    cox_results <- fit_weighted_cox(Y_combined, delta_combined, Z_combined,
                                    weights_combined, n_boot = n_boot,
                                    X_internal = X_internal,
                                    X_external = X_external,
                                    weight_type = weight_type,
                                    weight_params = weight_params)

    # Calculate performance metrics
    ess <- calculate_ess(weights_internal)
    max_weight_val <- max(weights_internal)

    # Assess balance
    # For under-specified, still check balance on all prognostic factors
    if (weighting_scenario == "correct") {
      balance_X <- X_internal[, 1:3]
      balance_means <- colMeans(X_external[, 1:3])
    } else if (weighting_scenario == "under") {
      # Check balance on all prognostic factors (including omitted X3)
      balance_X <- X_internal[, 1:3]
      balance_means <- colMeans(X_external[, 1:3])
    } else {
      balance_X <- X_internal
      balance_means <- colMeans(X_external)
    }

    max_std_diff <- calculate_balance(balance_X, balance_means, weights_internal)

    # Store results for this method
    results <- rbind(results, data.frame(
      n = n,
      conditional_hr = conditional_hr,
      true_marginal_hr = true_marginal_hr,
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
      max_weight = max_weight_val,
      max_std_diff = max_std_diff,
      n_boot_valid = cox_results$n_boot_valid,
      stringsAsFactors = FALSE
    ))
  }

  # ===========================================
  # STEP 6: Process IPTW ATUT Methods
  # ===========================================

  iptw_methods <- c("IPTW_OW", "IPTW_SW1", "IPTW_SWN", "IPTW_SWESS")

  for (method in iptw_methods) {
    # Extract normalization from method name
    norm_method <- substring(method, 6)  # Remove "IPTW_" prefix

    # Normalize IPTW weights for INTERNAL group
    weights_internal <- normalize_weights(iptw_weights, norm_method)

    # Combined weights: Internal = weighted, External = 1
    weights_combined <- c(weights_internal, rep(1, n))

    # Prepare parameters
    weight_params <- list(
      norm_method = norm_method,
      max_weight = max_weight,
      X_for_ps = X_for_ps
    )

    # Fit Cox model with bootstrap
    cox_results <- fit_weighted_cox(Y_combined, delta_combined, Z_combined,
                                    weights_combined, n_boot = n_boot,
                                    X_internal = X_internal,
                                    X_external = X_external,
                                    weight_type = "iptw",
                                    weight_params = weight_params)

    # Calculate metrics
    ess <- calculate_ess(weights_internal)
    max_weight_val <- max(weights_internal)

    # Balance assessment (same as above)
    if (weighting_scenario == "correct") {
      balance_X <- X_internal[, 1:3]
      balance_means <- colMeans(X_external[, 1:3])
    } else if (weighting_scenario == "under") {
      balance_X <- X_internal[, 1:3]
      balance_means <- colMeans(X_external[, 1:3])
    } else {
      balance_X <- X_internal
      balance_means <- colMeans(X_external)
    }

    max_std_diff <- calculate_balance(balance_X, balance_means, weights_internal)

    # Store IPTW results
    results <- rbind(results, data.frame(
      n = n,
      conditional_hr = conditional_hr,
      true_marginal_hr = true_marginal_hr,
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
      max_weight = max_weight_val,
      max_std_diff = max_std_diff,
      n_boot_valid = cox_results$n_boot_valid,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}

#' Run full simulation study
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
    conditional_hr = sim_params$conditional_hr,
    distribution = sim_params$distributions,
    median_surv = sim_params$median_surv,
    surv_48m = sim_params$surv_48m,
    dropout_rate = sim_params$dropout_rates,
    event_rate = sim_params$event_rates,
    ess_ratio = sim_params$ess_ratios,
    weighting_scenario = c("correct", "under", "over"),
    stringsAsFactors = FALSE
  )

  scenarios$d <- d_values[as.character(scenarios$ess_ratio)]
  scenarios$max_followup <- sim_params$max_followup

  # Initialize columns for beta and true_marginal_hr
  scenarios$beta1 <- NA
  scenarios$beta2 <- NA
  scenarios$beta3 <- NA
  scenarios$beta4 <- NA
  scenarios$gamma <- NA
  scenarios$true_marginal_hr <- NA

  # Calibrate beta, gamma, and compute true marginal HR
  cat("Calibrating parameters and computing true marginal HRs...\n")

  unique_combos <- unique(scenarios[, c("distribution", "median_surv", "surv_48m",
                                        "dropout_rate", "event_rate", "conditional_hr",
                                        "d")])

  for (i in 1:nrow(unique_combos)) {
    combo <- unique_combos[i, ]
    cat(sprintf("  Combo %d/%d: ", i, nrow(unique_combos)))

    base_params <- calibrate_base_params(combo$distribution, combo$median_surv,
                                         combo$surv_48m, target_time = sim_params$max_followup)

    beta <- calibrate_beta(combo$event_rate, distribution = combo$distribution,
                           base_params = base_params, dropout_rate = combo$dropout_rate,
                           max_followup = sim_params$max_followup)

    gamma <- log(combo$conditional_hr)

    # Compute true marginal ATUT HR
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

    cat(sprintf("Conditional HR = %.3f, Marginal HR = %.3f\n",
                combo$conditional_hr, true_marginal))

    # Find matching rows and update
    matching_rows <- which(
      scenarios$distribution == combo$distribution &
        scenarios$median_surv == combo$median_surv &
        scenarios$surv_48m == combo$surv_48m &
        scenarios$dropout_rate == combo$dropout_rate &
        scenarios$event_rate == combo$event_rate &
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

  # Export all necessary functions
  clusterExport(cl, c("generate_covariates", "calibrate_base_params",
                      "generate_survival_times", "generate_censoring",
                      "calculate_maic_weights", "calculate_mew_weights",
                      "calculate_iptw_atut_weights",
                      "normalize_weights", "calculate_ess", "fit_weighted_cox",
                      "calculate_balance", "run_single_simulation_with_params",
                      "compute_true_marginal_atut_hr"),
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

    # Reconstruct beta vector from individual columns
    current_scenario$beta <- c(
      current_scenario$beta1,
      current_scenario$beta2,
      current_scenario$beta3,
      current_scenario$beta4
    )

    # Create parameter list for each replication
    param_list <- replicate(sim_params$n_sim, current_scenario, simplify = FALSE)

    # Run parallel simulations
    scenario_results <- parLapply(cl, param_list, function(params) {
      tryCatch({
        run_single_simulation_with_params(params, sim_params)
      }, error = function(e) {
        cat("\nError in simulation:", e$message, "\n")
        return(NULL)
      })
    })

    # Filter out NULL results
    scenario_results <- scenario_results[!sapply(scenario_results, is.null)]

    if (length(scenario_results) > 0) {
      scenario_results_df <- do.call(rbind, scenario_results)
      scenario_results_df$iteration <- rep(1:length(scenario_results),
                                           each = 12)  # 12 methods
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

  # Combine all results
  final_results <- do.call(rbind, all_results)

  # Save final results
  if (save_intermediate) {
    filename_base <- sprintf("%s_n%d-%d_chr%.2f-%.2f_sim%d",
                             sim_params$results_prefix,
                             min(sim_params$n_sizes), max(sim_params$n_sizes),
                             min(sim_params$conditional_hr), max(sim_params$conditional_hr),
                             sim_params$n_sim)

    saveRDS(final_results, file.path(save_dir, paste0(filename_base, ".rds")))
    write.csv(final_results, file.path(save_dir, paste0(filename_base, ".csv")),
              row.names = FALSE)

    # Save scenario summary
    write.csv(scenarios, file.path(save_dir, "scenario_definitions.csv"),
              row.names = FALSE)

    cat(sprintf("\nResults saved to %s/\n", save_dir))
  }

  return(final_results)
}



# =============================================================================
# PERFORMANCE EVALUATION FUNCTIONS
# =============================================================================

#' Calculate performance metrics
#'
#' @description Computes bias, coverage, and efficiency metrics
#' Uses TRUE MARGINAL HR as the truth (not conditional HR)
calculate_performance_metrics <- function(results) {
  cat("Calculating performance metrics using true marginal ATUT HR...\n")

  # Calculate bias using true marginal HR
  results$bias <- results$hr_est - results$true_marginal_hr
  results$rel_bias <- (results$hr_est - results$true_marginal_hr) / results$true_marginal_hr
  results$abs_bias <- abs(results$bias)
  results$abs_rel_bias <- abs(results$rel_bias)

  # Coverage indicators using true marginal HR
  results$covered_robust <- (results$true_marginal_hr >= results$ci_robust_lower) &
    (results$true_marginal_hr <= results$ci_robust_upper)
  results$covered_boot_perc <- (results$true_marginal_hr >= results$ci_boot_perc_lower) &
    (results$true_marginal_hr <= results$ci_boot_perc_upper)
  results$covered_boot_norm <- (results$true_marginal_hr >= results$ci_boot_norm_lower) &
    (results$true_marginal_hr <= results$ci_boot_norm_upper)
  results$covered_boot_bca <- (results$true_marginal_hr >= results$ci_boot_bca_lower) &
    (results$true_marginal_hr <= results$ci_boot_bca_upper)

  # CI widths
  results$ci_width_robust <- results$ci_robust_upper - results$ci_robust_lower
  results$ci_width_boot_perc <- results$ci_boot_perc_upper - results$ci_boot_perc_lower
  results$ci_width_boot_norm <- results$ci_boot_norm_upper - results$ci_boot_norm_lower
  results$ci_width_boot_bca <- results$ci_boot_bca_upper - results$ci_boot_bca_lower

  # Summary by scenario and method
  summary_metrics <- aggregate(
    cbind(bias, rel_bias, abs_bias, abs_rel_bias,
          covered_robust, covered_boot_perc, covered_boot_norm, covered_boot_bca,
          ci_width_robust, ci_width_boot_perc, ci_width_boot_norm, ci_width_boot_bca,
          ess, max_weight, max_std_diff) ~
      n + conditional_hr + true_marginal_hr + distribution + weighting_scenario + method,
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

  # Add sample size information
  count_summary <- aggregate(
    hr_est ~ n + conditional_hr + true_marginal_hr + distribution + weighting_scenario + method,
    data = results,
    FUN = function(x) c(n_valid = sum(!is.na(x)), n_total = length(x))
  )

  summary_metrics <- merge(summary_metrics, count_summary,
                           by = c("n", "conditional_hr", "true_marginal_hr",
                                  "distribution", "weighting_scenario", "method"))

  # Print key findings
  cat("\n=== Performance Summary ===\n")
  cat("Bias (using true marginal HR):\n")
  bias_summary <- aggregate(bias ~ method, data = results,
                            FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                sd = sd(x, na.rm = TRUE)))
  print(bias_summary)

  cat("\nCoverage (95% CI):\n")
  coverage_summary <- aggregate(covered_robust ~ method, data = results,
                                FUN = function(x) mean(x, na.rm = TRUE))
  print(coverage_summary)

  cat("\nNon-collapsibility check:\n")
  nc_check <- unique(results[, c("conditional_hr", "true_marginal_hr")])
  nc_check$ratio <- nc_check$true_marginal_hr / nc_check$conditional_hr
  print(nc_check)

  return(summary_metrics)
}

#' Create comparison tables for MAIC vs IPTW
#'
#' @description Compares performance of MAIC/MEW vs IPTW ATUT methods
create_comparison_tables <- function(results) {
  cat("\nCreating MAIC/MEW vs IPTW ATUT comparison tables...\n")

  # Separate method types
  results$method_type <- ifelse(grepl("^IPTW", results$method), "IPTW_ATUT",
                                ifelse(grepl("^MEW", results$method), "MEW", "MAIC"))

  # Calculate metrics using true marginal HR
  results$bias <- results$hr_est - results$true_marginal_hr
  results$rel_bias <- (results$hr_est - results$true_marginal_hr) / results$true_marginal_hr
  results$mse <- (results$hr_est - results$true_marginal_hr)^2

  # Overall comparison
  overall_comparison <- aggregate(
    cbind(bias, rel_bias, mse, ess) ~ method_type + n + weighting_scenario,
    data = results,
    FUN = function(x) {
      c(mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE))
    }
  )

  # Best method by scenario
  best_by_scenario <- aggregate(
    mse ~ method + n + conditional_hr + weighting_scenario,
    data = results,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  best_method <- aggregate(
    mse ~ n + conditional_hr + weighting_scenario,
    data = best_by_scenario,
    FUN = function(x) {
      idx <- which.min(x)
      return(idx)
    }
  )

  # Relative efficiency (IPTW as reference)
  iptw_reference <- results[results$method == "IPTW_SW1", ]

  efficiency_comparison <- list()
  for (m in unique(results$method)) {
    if (m != "IPTW_SW1") {
      method_data <- results[results$method == m, ]
      merged <- merge(method_data, iptw_reference,
                      by = c("n", "conditional_hr", "distribution",
                             "weighting_scenario", "iteration"),
                      suffixes = c("", "_iptw"))

      merged$relative_mse <- merged$mse / merged$mse_iptw
      merged$relative_variance <- (merged$hr_est - mean(merged$hr_est, na.rm = TRUE))^2 /
        (merged$hr_est_iptw - mean(merged$hr_est_iptw, na.rm = TRUE))^2

      efficiency_comparison[[m]] <- data.frame(
        method = m,
        mean_relative_mse = mean(merged$relative_mse, na.rm = TRUE),
        median_relative_mse = median(merged$relative_mse, na.rm = TRUE),
        mean_relative_var = mean(merged$relative_variance, na.rm = TRUE)
      )
    }
  }

  efficiency_df <- do.call(rbind, efficiency_comparison)

  return(list(
    overall = overall_comparison,
    best_method = best_method,
    efficiency = efficiency_df
  ))
}

#' Create summary tables
#'
#' @description Creates formatted summary tables for reporting
create_summary_tables <- function(results, output_format = "latex") {
  cat("\nCreating summary tables...\n")

  # Main results table
  main_table <- aggregate(
    cbind(bias, rel_bias, covered_robust, ci_width_robust, ess) ~
      method + n + weighting_scenario,
    data = results,
    FUN = function(x) {
      c(mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE))
    }
  )

  # Format for output
  if (output_format == "latex") {
    # Format as LaTeX table
    cat("\n\\begin{table}[h]\n")
    cat("\\caption{Performance metrics by method and sample size}\n")
    cat("\\begin{tabular}{llrrrrrr}\n")
    cat("\\toprule\n")
    cat("Method & n & Bias & Rel. Bias & Coverage & CI Width & ESS \\\\\n")
    cat("\\midrule\n")

    for (i in 1:nrow(main_table)) {
      row <- main_table[i, ]
      cat(sprintf("%s & %d & %.3f (%.3f) & %.3f (%.3f) & %.3f & %.3f (%.3f) & %.1f (%.1f) \\\\\n",
                  row$method, row$n,
                  row$bias[1], row$bias[2],
                  row$rel_bias[1], row$rel_bias[2],
                  row$covered_robust[1],
                  row$ci_width_robust[1], row$ci_width_robust[2],
                  row$ess[1], row$ess[2]))
    }

    cat("\\bottomrule\n")
    cat("\\end{tabular}\n")
    cat("\\end{table}\n")
  }

  return(main_table)
}

# =============================================================================
# ANALYSIS HELPER FUNCTIONS
# =============================================================================

#' Check weight distributions across methods
#'
#' @description Analyzes weight properties (ESS, max weight) by method and scenario
#' @param results Simulation results dataframe
#' @return Summary dataframe with weight statistics
check_weight_distributions <- function(results) {
  cat("Analyzing weight distributions...\n")

  weight_summary <- aggregate(
    cbind(ess, max_weight) ~ method + weighting_scenario,
    data = results,
    FUN = function(x) {
      if (all(is.na(x))) {
        return(c(mean = NA, sd = NA, min = NA, max = NA, median = NA))
      }
      c(mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        min = min(x, na.rm = TRUE),
        max = max(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE))
    }
  )

  # Add method type classification
  weight_summary$method_type <- ifelse(grepl("^IPTW", weight_summary$method), "IPTW_ATUT",
                                       ifelse(grepl("^MEW", weight_summary$method), "MEW", "MAIC"))

  return(weight_summary)
}

#' Validate bootstrap results
#'
#' @description Checks bootstrap convergence and validity rates
validate_bootstrap <- function(results, expected_n_boot = NULL) {
  cat("Validating bootstrap results...\n")

  if (is.null(expected_n_boot)) {
    expected_n_boot <- max(results$n_boot_valid, na.rm = TRUE)
    cat(sprintf("  Inferred n_boot = %d\n", expected_n_boot))
  }

  boot_validity <- aggregate(
    n_boot_valid ~ method + n + weighting_scenario,
    data = results,
    FUN = function(x) {
      valid_prop <- mean(x / expected_n_boot, na.rm = TRUE)
      c(mean_valid = mean(x, na.rm = TRUE),
        prop_valid = valid_prop,
        min_valid = min(x, na.rm = TRUE),
        max_valid = max(x, na.rm = TRUE),
        n_failed = sum(x < expected_n_boot * 0.5))
    }
  )

  boot_validity$flag_low_validity <- boot_validity$n_boot_valid[, "prop_valid"] < 0.8

  return(boot_validity)
}

#' Compare MAIC/MEW vs IPTW ATUT performance
#'
#' @description Direct comparison of method types
compare_methods <- function(results) {
  cat("Comparing MAIC/MEW vs IPTW ATUT methods...\n")

  # Separate method types
  maic_methods <- results[!grepl("^IPTW", results$method), ]
  iptw_methods <- results[grepl("^IPTW", results$method), ]

  # Calculate bias using true marginal HR
  maic_methods$bias <- maic_methods$hr_est - maic_methods$true_marginal_hr
  iptw_methods$bias <- iptw_methods$hr_est - iptw_methods$true_marginal_hr

  # Calculate relative bias
  maic_methods$rel_bias <- maic_methods$bias / maic_methods$true_marginal_hr
  iptw_methods$rel_bias <- iptw_methods$bias / iptw_methods$true_marginal_hr

  # Summary statistics
  comparison <- data.frame(
    Method_Type = c("MAIC/MEW", "IPTW_ATUT"),
    N_Obs = c(nrow(maic_methods), nrow(iptw_methods)),
    Mean_HR = c(mean(maic_methods$hr_est, na.rm = TRUE),
                mean(iptw_methods$hr_est, na.rm = TRUE)),
    Mean_Bias = c(mean(maic_methods$bias, na.rm = TRUE),
                  mean(iptw_methods$bias, na.rm = TRUE)),
    Mean_Abs_Bias = c(mean(abs(maic_methods$bias), na.rm = TRUE),
                      mean(abs(iptw_methods$bias), na.rm = TRUE)),
    Mean_Rel_Bias = c(mean(maic_methods$rel_bias, na.rm = TRUE),
                      mean(iptw_methods$rel_bias, na.rm = TRUE)),
    SD_HR = c(sd(maic_methods$hr_est, na.rm = TRUE),
              sd(iptw_methods$hr_est, na.rm = TRUE)),
    Mean_ESS = c(mean(maic_methods$ess, na.rm = TRUE),
                 mean(iptw_methods$ess, na.rm = TRUE)),
    Mean_Max_Weight = c(mean(maic_methods$max_weight, na.rm = TRUE),
                        mean(iptw_methods$max_weight, na.rm = TRUE)),
    Mean_Balance = c(mean(maic_methods$max_std_diff, na.rm = TRUE),
                     mean(iptw_methods$max_std_diff, na.rm = TRUE)),
    Mean_Bootstrap_Valid = c(mean(maic_methods$n_boot_valid, na.rm = TRUE),
                             mean(iptw_methods$n_boot_valid, na.rm = TRUE)),
    Prop_NA = c(mean(is.na(maic_methods$hr_est)),
                mean(is.na(iptw_methods$hr_est)))
  )

  cat("\n=== Key Comparison: MAIC/MEW vs IPTW ATUT ===\n")
  print(comparison)

  return(comparison)
}

#' Check simulation convergence
#'
#' @description Assesses Monte Carlo error and convergence
summarize_convergence <- function(results, n_sim = NULL) {
  cat("Checking simulation convergence...\n")

  # Count iterations per scenario
  scenario_counts <- aggregate(
    iteration ~ n + conditional_hr + distribution + weighting_scenario + method,
    data = results,
    FUN = function(x) length(unique(x))
  )

  names(scenario_counts)[names(scenario_counts) == "iteration"] <- "n_iterations"

  if (!is.null(n_sim)) {
    scenario_counts$completion_rate <- scenario_counts$n_iterations / n_sim
    scenario_counts$flag_incomplete <- scenario_counts$completion_rate < 0.9
  }

  # Calculate Monte Carlo standard errors using true marginal HR
  results$squared_error <- (results$hr_est - results$true_marginal_hr)^2

  mc_errors <- aggregate(
    cbind(hr_est, squared_error) ~ n + conditional_hr + distribution +
      weighting_scenario + method,
    data = results,
    FUN = function(x) {
      n_valid <- sum(!is.na(x))
      if (n_valid > 1) {
        se <- sd(x, na.rm = TRUE) / sqrt(n_valid)
        cv <- sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
        return(c(n_valid = n_valid,
                 mean = mean(x, na.rm = TRUE),
                 sd = sd(x, na.rm = TRUE),
                 se = se,
                 cv = cv))
      } else {
        return(c(n_valid = n_valid, mean = NA, sd = NA, se = NA, cv = NA))
      }
    }
  )

  # Merge results
  convergence_summary <- merge(scenario_counts, mc_errors,
                               by = c("n", "conditional_hr", "distribution",
                                      "weighting_scenario", "method"))

  # Flag high MC error
  convergence_summary$flag_high_mc_error <-
    convergence_summary$hr_est[, "se"] > 0.05 * abs(convergence_summary$hr_est[, "mean"])

  # Summary statistics
  cat(sprintf("\nTotal scenarios: %d\n", nrow(convergence_summary)))

  if (!is.null(n_sim)) {
    cat(sprintf("Scenarios with <90%% completion: %d\n",
                sum(convergence_summary$flag_incomplete, na.rm = TRUE)))
  }

  cat(sprintf("Scenarios with high MC error: %d\n",
              sum(convergence_summary$flag_high_mc_error, na.rm = TRUE)))

  return(convergence_summary)
}

#' Main function with integrated analysis
#'
#' @param sim_params List containing all simulation parameters
#' @param run_analysis Whether to run analysis helpers after simulation
#' @return List with results, metrics, and analysis summaries
main <- function(sim_params = NULL, run_analysis = TRUE) {
  if (is.null(sim_params)) {
    stop("sim_params must be provided")
  }

  # Set defaults
  if (is.null(sim_params$timestamp_files)) sim_params$timestamp_files <- TRUE
  if (is.null(sim_params$results_prefix)) sim_params$results_prefix <- "maic_iptw_atut_sim"
  if (is.null(sim_params$n_boot)) sim_params$n_boot <- 500
  if (is.null(sim_params$max_weight)) sim_params$max_weight <- 100
  if (is.null(sim_params$n_mc_truth)) sim_params$n_mc_truth <- 10000

  start_time <- Sys.time()

  cat("=============================================================================\n")
  cat("MAIC Simulation Study with IPTW ATUT Comparison - Version 39\n")
  cat("TARGET ESTIMAND: ATUT (Average Treatment Effect in the Untreated)\n")
  cat("=============================================================================\n\n")

  # Run simulation
  cat("Starting simulation study...\n")
  results <- run_simulation_study(sim_params)

  # Initialize analysis results list
  analysis_results <- list()

  if (run_analysis && !is.null(results) && nrow(results) > 0) {
    cat("\n=== Running Post-Simulation Analyses ===\n")

    # 1. Calculate performance metrics
    cat("\n1. Performance Metrics (using true marginal ATUT HR)\n")
    analysis_results$performance_metrics <- tryCatch({
      calculate_performance_metrics(results)
    }, error = function(e) {
      cat("  Error in performance metrics:", e$message, "\n")
      NULL
    })

    # 2. Check weight distributions
    cat("\n2. Weight Distribution Analysis\n")
    analysis_results$weight_distributions <- tryCatch({
      check_weight_distributions(results)
    }, error = function(e) {
      cat("  Error in weight analysis:", e$message, "\n")
      NULL
    })

    # 3. Validate bootstrap results
    cat("\n3. Bootstrap Validation\n")
    analysis_results$bootstrap_validation <- tryCatch({
      validate_bootstrap(results, expected_n_boot = sim_params$n_boot)
    }, error = function(e) {
      cat("  Error in bootstrap validation:", e$message, "\n")
      NULL
    })

    # 4. Compare methods
    cat("\n4. Method Comparison (MAIC/MEW vs IPTW ATUT)\n")
    analysis_results$method_comparison <- tryCatch({
      compare_methods(results)
    }, error = function(e) {
      cat("  Error in method comparison:", e$message, "\n")
      NULL
    })

    # 5. Create comparison tables
    cat("\n5. Detailed Comparison Tables\n")
    analysis_results$comparison_tables <- tryCatch({
      create_comparison_tables(results)
    }, error = function(e) {
      cat("  Error in comparison tables:", e$message, "\n")
      NULL
    })

    # 6. Check convergence
    cat("\n6. Convergence Analysis\n")
    analysis_results$convergence <- tryCatch({
      summarize_convergence(results, n_sim = sim_params$n_sim)
    }, error = function(e) {
      cat("  Error in convergence analysis:", e$message, "\n")
      NULL
    })

    # Print key summaries
    if (!is.null(analysis_results$method_comparison)) {
      cat("\n=== Method Comparison Summary ===\n")
      print(analysis_results$method_comparison)
    }

    if (!is.null(analysis_results$comparison_tables)) {
      cat("\n=== Relative Efficiency (vs IPTW_SW1) ===\n")
      print(analysis_results$comparison_tables$efficiency[1:10, ])  # Show top 10
    }

    # Check for non-collapsibility
    cat("\n=== Non-collapsibility Assessment ===\n")
    nc_data <- unique(results[, c("conditional_hr", "true_marginal_hr")])
    nc_data$ratio <- nc_data$true_marginal_hr / nc_data$conditional_hr
    nc_data$percent_diff <- (nc_data$true_marginal_hr - nc_data$conditional_hr) /
      nc_data$conditional_hr * 100
    print(nc_data)

    if (!is.null(analysis_results$bootstrap_validation)) {
      boot_summary <- analysis_results$bootstrap_validation
      low_validity <- sum(boot_summary$flag_low_validity, na.rm = TRUE)
      if (low_validity > 0) {
        cat(sprintf("\nWarning: %d scenarios have <80%% bootstrap validity\n", low_validity))
      }
    }

    if (!is.null(analysis_results$convergence)) {
      conv_summary <- analysis_results$convergence
      high_error <- sum(conv_summary$flag_high_mc_error, na.rm = TRUE)
      if (high_error > 0) {
        cat(sprintf("\nWarning: %d scenarios have high Monte Carlo error\n", high_error))
      }
    }
  }

  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "hours")

  cat("\n=============================================================================\n")
  cat(sprintf("Simulation completed in %.2f hours\n", runtime))
  cat(sprintf("Total scenarios: %d\n", length(unique(paste(results$n, results$conditional_hr,
                                                           results$distribution,
                                                           results$weighting_scenario)))))
  cat(sprintf("Total observations: %d\n", nrow(results)))
  cat(sprintf("Methods evaluated: %d\n", length(unique(results$method))))
  cat("=============================================================================\n")

  # Save analysis results if requested
  if (sim_params$timestamp_files) {
    save_dir <- paste0("simulation_results_", format(start_time, "%Y%m%d_%H%M%S"))
  } else {
    save_dir <- "simulation_results"
  }

  if (dir.exists(save_dir) && !is.null(analysis_results)) {
    saveRDS(analysis_results, file.path(save_dir, "analysis_summaries.rds"))
    cat(sprintf("\nAnalysis summaries saved to %s/analysis_summaries.rds\n", save_dir))
  }

  return(list(
    results = results,
    analysis = analysis_results,
    runtime = runtime,
    sim_params = sim_params
  ))
}
