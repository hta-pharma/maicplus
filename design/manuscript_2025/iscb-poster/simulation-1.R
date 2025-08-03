# =============================================================================
# MAIC Simulation Study for Time-to-Event Outcomes (Version A)
#
# This code implements the simulation study described in the manuscript:
# "Comparative Evaluation of Weight Normalization Methods in Matching-Adjusted
# Indirect Comparison for Time-to-Event Outcomes"
#
# Developer: Gregory etc.
# Date: Aug, 2025
# R Version: 4.5.1
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
#' @return Calibrated d value
calibrate_d_for_ess <- function(target_ess_ratio, n_calib = 10000, tol = 0.01) {
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

    # Compute ESS
    ess <- sum(weights)^2 / sum(weights^2)

    # Return ESS ratio
    return(ess / (2 * n_calib))
  }

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

  return((d_lower + d_upper) / 2)
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
