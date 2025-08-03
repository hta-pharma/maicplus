# =============================================================================
# MAIC Simulation Study for Time-to-Event Outcomes (Version B)
#
# This code implements the simulation study described in the manuscript:
# "Comparative Evaluation of Weight Normalization Methods in Matching-Adjusted
# Indirect Comparison for Time-to-Event Outcomes"
#
# Developer: Gregory etc.
# Date: Aug, 2025
# R Version: 4.5.1
# =============================================================================

# Load minimal required packages
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Please install the 'survival' package to proceed.")
}
library(survival)
library(parallel)

# ------------------------
# Utility functions
# ------------------------

# Effective Sample Size (ESS) computation
compute_ESS <- function(w) {
  # w: vector of (possibly unnormalized) weights
  num <- sum(w)^2
  den <- sum(w^2)
  if (den <= 0) return(0)
  return(num / den)
}

# Solve exponential tilting (MAIC) weights to match target means on selected covariates
# Implements Newton-Raphson for solving: weighted_means(theta) = target_means
solve_tilting <- function(X, target_means, max_iter = 100, tol = 1e-8) {
  # X: matrix (n x p) of covariates to match (each column j)
  # target_means: length-p vector of desired weighted means
  n <- nrow(X)
  p <- ncol(X)
  theta <- rep(0, p)  # initial guess
  for (iter in seq_len(max_iter)) {
    # Compute weights unnormalized: w_i = exp(theta^T x_i)
    eta <- as.vector(X %*% theta)
    w_unnorm <- exp(eta)
    sum_w <- sum(w_unnorm)
    # Weighted means
    weighted_means <- colSums(w_unnorm * X) / sum_w
    # Residual
    g <- weighted_means - target_means
    if (max(abs(g)) < tol) break
    # Compute covariance matrix under current weights: Cov_w(X)
    # Cov_w(X) = E_w[XX^T] - (E_w[X]) (E_w[X])^T
    EwX <- weighted_means
    EwXX <- t(X) %*% (w_unnorm * X) / sum_w  # p x p matrix
    cov_w <- EwXX - tcrossprod(EwX)
    # Jacobian of g is cov_w
    # Update: theta_new = theta - inv(cov_w) %*% g
    # Regularize if cov_w is near-singular
    # Use tryCatch to avoid failures
    delta_theta <- tryCatch({
      solve(cov_w, g)
    }, error = function(e) {
      # fallback: small ridge
      solve(cov_w + diag(1e-6, p), g)
    })
    theta <- theta - delta_theta
    # Damping to keep stable (optional)
    if (max(abs(delta_theta)) < tol) break
  }
  # Final weights normalized to sum to n (matching original scaling used in manuscript)
  w_raw <- exp(as.vector(X %*% theta))
  w <- w_raw * (n / sum(w_raw))  # sum-to-N normalization
  return(list(weights = w, theta = theta, converged = (max(abs(g)) < tol)))
}

# MAIC weight computation wrapper: optionally include nuisance covariate
maic_weights <- function(X_full, target_means, prognostic_idx = 1:3, include_nuisance = FALSE) {
  # X_full: full covariate matrix (n x 4), columns correspond to X1..X4
  # target_means: length-3 or length-4 depending on include_nuisance
  if (include_nuisance) {
    X <- X_full  # match on all provided (including X4)
  } else {
    X <- X_full[, prognostic_idx, drop = FALSE]  # only prognostic
  }
  res <- solve_tilting(X, target_means = target_means)
  return(res$weights)
}

# Standardized difference between two samples for covariate j
std_diff <- function(x_weighted, x_ref, w = NULL) {
  # x_weighted: vector of covariate in weighted sample (e.g., treated after weighting)
  # x_ref: vector of covariate in reference (control) sample
  # w: weights for x_weighted (if NULL assumes equal weights)
  if (is.null(w)) {
    m1 <- mean(x_weighted)
    s1 <- var(x_weighted)
  } else {
    m1 <- sum(w * x_weighted) / sum(w)
    # Weighted variance (biased version, approximated)
    mu <- m1
    s1 <- sum(w * (x_weighted - mu)^2) / sum(w)
  }
  m0 <- mean(x_ref)
  s0 <- var(x_ref)
  # Pooled sd (using control and weighted sample) or use control-only SD
  pooled_sd <- sqrt((s0 + s1)/2)
  if (pooled_sd == 0) return(0)
  return((m1 - m0) / pooled_sd)
}

# Generate baseline covariates for one population given mean vector
generate_baseline <- function(n, mu_vec) {
  # mu_vec: length-4 vector of means; for continuous variable X1, this is its mean; for binaries X2-X4, success probs
  X1 <- rnorm(n, mean = mu_vec[1], sd = 0.5)
  X2 <- rbinom(n, size = 1, prob = min(max(mu_vec[2], 0), 1))
  X3 <- rbinom(n, size = 1, prob = min(max(mu_vec[3], 0), 1))
  X4 <- rbinom(n, size = 1, prob = min(max(mu_vec[4], 0), 1))
  X <- cbind(X1, X2, X3, X4)
  colnames(X) <- paste0("X", 1:4)
  return(X)
}

# Calibration of delta (scaling on mu0) to achieve a target overlap ESS/N = p
calibrate_delta <- function(target_p, mu0 = c(0.5, 0.5, 0.5, 0.5), prognostic_idx = 1:3,
                            Nlarge = 20000, tol = 1e-3, max_iter = 30) {
  # We reinterpret: external population has mean mu1 = delta * mu0, and we weight external to internal (reference mu0)
  f <- function(delta) {
    mu1 <- delta * mu0
    # Draw large external sample with means mu1, and internal reference with mu0
    X_ext <- generate_baseline(Nlarge, mu1)
    X_ref <- generate_baseline(Nlarge, mu0)
    # Compute MAIC weights: re-weight external to match prognostic means of internal
    target_means <- colMeans(X_ref[, prognostic_idx, drop = FALSE])
    w_ext <- maic_weights(X_ext, target_means = target_means, prognostic_idx = prognostic_idx, include_nuisance = FALSE)
    ess <- compute_ESS(w_ext)
    ratio <- ess / Nlarge
    return(ratio - target_p)
  }
  # Bracket search: delta in (eps, 2]; ensure sign change
  lower <- 0.01
  upper <- 2
  f_lower <- f(lower)
  f_upper <- f(upper)
  if (f_lower * f_upper > 0) {
    warning("Unable to bracket root in calibrate_delta; returning delta=1 (no scaling).")
    return(1)
  }
  # Use uniroot to solve for delta
  sol <- uniroot(f, lower = lower, upper = upper, tol = tol, maxiter = max_iter)
  return(sol$root)
}

# Compute baseline Weibull parameters (shape k, scale lambda) from median and survival at time t0
weibull_params_from_summary <- function(median, surv_at_48, t0 = 48) {
  # Solve for shape k from:
  # surv_at_48 = exp( - (t0 / lambda)^k )
  # median = lambda * (log 2)^{1/k}
  # => k = log( -log(surv_at_48)/log(2) ) / log(t0 / median)
  numerator <- log(-log(surv_at_48) / log(2))
  denominator <- log(t0 / median)
  k <- numerator / denominator
  lambda <- median / (log(2))^(1 / k)
  return(list(shape = k, scale = lambda))
}

# Compute log-normal parameters (mu, sigma) from median and survival at t0=48
lognormal_params_from_summary <- function(median, surv_at_48, t0 = 48) {
  # median => mu = log(median)
  mu <- log(median)
  # S(t0) = 1 - Phi( (log t0 - mu)/sigma ) = surv_at_48
  # => Phi( (log t0 - mu)/sigma ) = 1 - surv_at_48
  z <- qnorm(1 - surv_at_48)
  sigma <- (log(t0) - mu) / z
  return(list(mu = mu, sigma = sigma))
}

# Generate event times under proportional hazards for Weibull or log-normal baseline
generate_event_times <- function(X, Z, beta, gamma, dist = c("weibull", "lognormal"),
                                 median_surv = 24, surv48 = 0.05) {
  dist <- match.arg(dist)
  n <- nrow(X)
  # Linear predictor: X %*% beta + gamma * Z
  LP <- as.vector(X[, 1:length(beta), drop = FALSE] %*% beta + gamma * Z)
  if (dist == "weibull") {
    # Get baseline Weibull parameters
    wp <- weibull_params_from_summary(median_surv, surv48)
    k <- wp$shape
    lambda <- wp$scale
    # Generate uniform variates
    u <- runif(n)
    # T = lambda * ( -log(u) / exp(LP) )^{1/k}
    Ti <- lambda * ( -log(u) / exp(LP) )^(1 / k)
  } else if (dist == "lognormal") {
    lp <- lognormal_params_from_summary(median_surv, surv48)
    mu0 <- lp$mu
    sigma <- lp$sigma
    # Adjusted survival: use S(t) = S0(t)^{exp(LP)} => invert via
    # draw u ~ Unif(0,1), then baseline survival target v = u^{1 / exp(LP)}
    u <- runif(n)
    v <- u^(1 / exp(LP))
    # For log-normal, S0(t) = 1 - Phi((log t - mu)/sigma) = v
    # => Phi((log t - mu)/sigma) = 1 - v => log t = mu + sigma * qnorm(1 - v)
    logTi <- mu0 + sigma * qnorm(1 - v)
    Ti <- exp(logTi)
  } else {
    stop("Unsupported distribution.")
  }
  return(Ti)
}

# Simulate dropout times and apply censoring
apply_dropout_censoring <- function(T, P_drop) {
  n <- length(T)
  # total dropout fraction P_drop in {0.2, 0.4}
  # Early dropouts: 60% of P_drop in [0,6), rest 40% in [6,48]
  z_early <- rbinom(1, n, prob = 0.6 * P_drop)
  z_late <- rbinom(1, n, prob = 0.4 * P_drop)
  # Determine indices
  idx <- seq_len(n)
  drop_idx <- sample(idx, size = z_early + z_late, replace = FALSE)
  I_zeta06 <- drop_idx[seq_len(z_early)]
  I_zeta6p <- drop_idx[(z_early + 1):(z_early + z_late)]
  # Initialize censoring times
  C <- rep(48, n)  # administrative default
  # Early dropouts
  for (i in I_zeta06) {
    g1 <- min(T[i], 6)
    if (g1 > 0) {
      C[i] <- runif(1, 0, g1)
    } else {
      C[i] <- 0
    }
  }
  # Late dropouts
  for (i in I_zeta6p) {
    g2 <- min(T[i], 48)
    if (T[i] < 6) {
      # event occurred before 6: leave as administrative (48)
      C[i] <- 48
    } else {
      if (g2 > 6) {
        C[i] <- runif(1, 6, g2)
      } else {
        C[i] <- 48
      }
    }
  }
  # Dropout indicator zeta (1 if dropout happened before administrative)
  zeta <- as.integer(C < 48)
  Y <- pmin(T, C)
  delta <- as.integer(T <= C)
  return(list(Y = Y, delta = delta, censor_time = C, zeta = zeta))
}

# Estimate ATUT HR (treated reweighted to control) via weighted Cox model
estimate_ATUT_HR <- function(data, prognostic_idx = 1:3, include_nuisance = FALSE) {
  # data: data.frame containing columns Y, delta, Z, X1..X4
  # reference group is control Z=0; treated Z=1 will be reweighted to match controls on prognostic covariates
  X_full <- as.matrix(data[, paste0("X", 1:4)])
  control <- data[data$Z == 0, ]
  treated <- data[data$Z == 1, ]
  # Compute target means on prognostic covariates from control
  target_means <- colMeans(control[, paste0("X", prognostic_idx), drop = FALSE])
  # Compute MAIC weights for treated
  w_treated <- maic_weights(as.matrix(treated[, paste0("X", 1:4)]), target_means = target_means,
                            prognostic_idx = prognostic_idx, include_nuisance = include_nuisance)
  # Build full weight vector: controls weight=1, treated weights as computed
  w_full <- numeric(nrow(data))
  w_full[data$Z == 0] <- 1
  w_full[data$Z == 1] <- w_treated
  # Fit weighted Cox model with Z as covariate
  cox_mod <- tryCatch({
    coxph(Surv(Y, delta) ~ Z, data = data, weights = w_full, ties = "efron", control = coxph.control(iter.max = 50))
  }, error = function(e) {
    return(NULL)
  })
  if (is.null(cox_mod)) {
    return(list(hr = NA, ci = c(NA, NA)))
  }
  coef_est <- coef(cox_mod)["Z"]
  se <- sqrt(vcov(cox_mod)["Z", "Z"])
  hr <- exp(coef_est)
  ci_lower <- exp(coef_est - 1.96 * se)
  ci_upper <- exp(coef_est + 1.96 * se)
  return(list(hr = hr, ci = c(ci_lower, ci_upper), model = cox_mod, weights_treated = w_treated))
}

# Wrapper for calibration of gamma to hit target ATUT HR via Monte Carlo
calibrate_gamma <- function(target_hr, beta, dist = "weibull", median_surv = 24, surv48 = 0.05,
                            P_drop = 0.2, target_overlap_delta = 1, # baseline means scaling is assumed fixed
                            prognostic_idx = 1:3, Ncalib = 5000, tol = 1e-3, max_iter = 25, seed = 123) {
  # target_hr: desired marginal ATUT HR
  # beta: vector of length equal to number of prognostic covariates
  # target_overlap_delta: used if baseline mean vector is scaled; for simplicity assume control mean
  set.seed(seed)
  f <- function(gamma) {
    # Simulate large trial to estimate ATUT HR with given gamma
    n_each <- Ncalib
    # Generate baseline covariates for control (Z=0) and treated (Z=1)
    mu0 <- c(0.5, 0.5, 0.5, 0.5) * target_overlap_delta
    X_control <- generate_baseline(n_each, mu0)
    X_treated <- generate_baseline(n_each, mu0)
    # Assign Z
    Z_control <- rep(0, n_each)
    Z_treated <- rep(1, n_each)
    # Combine
    X_all <- rbind(X_control, X_treated)
    Z_all <- c(Z_control, Z_treated)
    # Generate event times
    T_control <- generate_event_times(X_control, Z_control, beta = beta, gamma = 0,
                                      dist = dist, median_surv = median_surv, surv48 = surv48)
    T_treated <- generate_event_times(X_treated, Z_treated, beta = beta, gamma = gamma,
                                      dist = dist, median_surv = median_surv, surv48 = surv48)
    # Apply dropout separately for each arm
    drop_control <- apply_dropout_censoring(T_control, P_drop)
    drop_treated <- apply_dropout_censoring(T_treated, P_drop)
    # Construct data
    data <- data.frame(
      Y = c(drop_control$Y, drop_treated$Y),
      delta = c(drop_control$delta, drop_treated$delta),
      Z = c(Z_control, Z_treated),
      X1 = c(X_control[, 1], X_treated[, 1]),
      X2 = c(X_control[, 2], X_treated[, 2]),
      X3 = c(X_control[, 3], X_treated[, 3]),
      X4 = c(X_control[, 4], X_treated[, 4])
    )
    est <- estimate_ATUT_HR(data, prognostic_idx = prognostic_idx, include_nuisance = FALSE)
    if (is.na(est$hr)) return(1e3)  # large penalty
    return(est$hr - target_hr)
  }
  # Bracket gamma in log-scale: exp(gamma) roughly between 0.1 and 2
  lower <- log(0.1)
  upper <- log(2)
  f_lower <- f(lower)
  f_upper <- f(upper)
  if (f_lower * f_upper > 0) {
    warning("Failed to bracket gamma for calibration; returning gamma = log(target_hr) as initial guess.")
    return(log(target_hr))
  }
  sol <- uniroot(f, lower = lower, upper = upper, tol = tol, maxiter = max_iter)
  return(sol$root)
}

# ------------------------
# Scenario runner
# ------------------------

# Single replicate of a scenario
run_single_rep <- function(n, true_hr, dist, median_surv, surv48, P_drop, beta, gamma,
                           delta_baseline, prognostic_spec = c("correct", "under", "over")) {
  # Generate baseline covariates: control and treated (same baseline for this example)
  mu0 <- c(0.5, 0.5, 0.5, 0.5)
  mu1 <- delta_baseline * mu0  # scaled mean for treated/internal depending on interpretation
  # For this replicate: use mu0 for control (Z=0) and mu1 for treated (Z=1) to reflect baseline differences
  X_control <- generate_baseline(n, mu0)
  X_treated <- generate_baseline(n, mu1)
  Z_control <- rep(0, n)
  Z_treated <- rep(1, n)
  # Generate event times
  T_control <- generate_event_times(X_control, Z_control, beta = beta, gamma = 0,
                                    dist = dist, median_surv = median_surv, surv48 = surv48)
  T_treated <- generate_event_times(X_treated, Z_treated, beta = beta, gamma = gamma,
                                    dist = dist, median_surv = median_surv, surv48 = surv48)
  # Apply dropout
  drop_control <- apply_dropout_censoring(T_control, P_drop)
  drop_treated <- apply_dropout_censoring(T_treated, P_drop)
  # Combine
  data <- data.frame(
    Y = c(drop_control$Y, drop_treated$Y),
    delta = c(drop_control$delta, drop_treated$delta),
    Z = c(Z_control, Z_treated),
    X1 = c(X_control[, 1], X_treated[, 1]),
    X2 = c(X_control[, 2], X_treated[, 2]),
    X3 = c(X_control[, 3], X_treated[, 3]),
    X4 = c(X_control[, 4], X_treated[, 4])
  )
  # Depending on prognostic_spec, decide which covariates to include in weighting
  res_list <- list()
  specs <- list(
    correct = list(include_nuisance = TRUE, prognostic_idx = 1:3),
    under = list(include_nuisance = FALSE, prognostic_idx = c(1, 2)), # omit X3 (key prognostic)
    over = list(include_nuisance = TRUE, prognostic_idx = 1:3) # includes nuisance but weights based on full or just prognostic? interpret as including nuisance in model
  )
  spec <- specs[[prognostic_spec]]
  include_nuisance <- (prognostic_spec == "correct" || prognostic_spec == "over")
  # Estimate ATUT HR under different specifications
  est <- estimate_ATUT_HR(data,
                          prognostic_idx = spec$prognostic_idx,
                          include_nuisance = include_nuisance)
  # Extract weight vector for treated if available
  # Compute balance metrics: standardized differences on X1..X4 between weighted treated and control
  # Recompute weights explicitly for treated (for diagnostics)
  control_data <- data[data$Z == 0, ]
  treated_data <- data[data$Z == 1, ]
  target_means <- colMeans(control_data[, paste0("X", spec$prognostic_idx), drop = FALSE])
  w_treated <- maic_weights(as.matrix(treated_data[, paste0("X", 1:4)]),
                            target_means = target_means,
                            prognostic_idx = spec$prognostic_idx,
                            include_nuisance = include_nuisance)
  # Weight normalizations (base = w_treated)
  OW <- w_treated
  SW1 <- OW / sum(OW)
  SWN <- OW * n / sum(OW)
  ESS_val <- compute_ESS(OW)
  SWESS <- OW * ESS_val / sum(OW)
  # Placeholder for MEW variants: here we fallback to base OW
  MEW <- OW
  MEW1 <- MEW / sum(MEW)
  MEWN <- MEW * n / sum(MEW)
  MEWESS <- MEW * ESS_val / sum(MEW)
  weight_matrix <- cbind(OW = OW, SW1 = SW1, SWN = SWN, SWESS = SWESS,
                         MEW = MEW, MEW1 = MEW1, MEWN = MEWN, MEWESS = MEWESS)
  # For each weight variant, compute HR via weighted Cox (treated reweighted; control weight=1)
  hr_estimates <- lapply(colnames(weight_matrix), function(wname) {
    w_full <- numeric(nrow(data))
    w_full[data$Z == 0] <- 1
    w_full[data$Z == 1] <- weight_matrix[, wname]
    cox_mod <- tryCatch({
      coxph(Surv(Y, delta) ~ Z, data = data, weights = w_full, ties = "efron",
            control = coxph.control(iter.max = 50))
    }, error = function(e) NULL)
    if (is.null(cox_mod)) return(list(hr = NA, ci = c(NA, NA)))
    coef_est <- coef(cox_mod)["Z"]
    se <- sqrt(vcov(cox_mod)["Z", "Z"])
    hr <- exp(coef_est)
    ci <- c(exp(coef_est - 1.96 * se), exp(coef_est + 1.96 * se))
    return(list(hr = hr, ci = ci))
  })
  names(hr_estimates) <- colnames(weight_matrix)
  # Compute bias, coverage for each variant
  results <- data.frame(
    method = colnames(weight_matrix),
    hr = sapply(hr_estimates, function(x) x$hr),
    ci_lower = sapply(hr_estimates, function(x) x$ci[1]),
    ci_upper = sapply(hr_estimates, function(x) x$ci[2]),
    bias = sapply(hr_estimates, function(x) x$hr - true_hr),
    coverage = sapply(hr_estimates, function(x) {
      ci <- x$ci
      if (is.na(ci[1]) || is.na(ci[2])) return(0)
      return(as.integer(ci[1] <= true_hr & true_hr <= ci[2]))
    }),
    ESS_treated = rep(ESS_val, length(colnames(weight_matrix))),
    stringsAsFactors = FALSE
  )
  # Balance metrics (standardized differences) after weighting treated toward control
  balance <- sapply(1:4, function(j) {
    x_control_j <- control_data[[paste0("X", j)]]
    x_treated_j <- treated_data[[paste0("X", j)]]
    std_diff(x_weighted = x_treated_j, x_ref = x_control_j, w = OW)
  })
  names(balance) <- paste0("std_diff_X", 1:4)
  results$std_diff_X1 <- balance["std_diff_X1"]
  results$std_diff_X2 <- balance["std_diff_X2"]
  results$std_diff_X3 <- balance["std_diff_X3"]
  results$std_diff_X4 <- balance["std_diff_X4"]
  # Return per-replicate metrics
  return(results)
}

# High-level scenario launcher (with parallelization)
run_scenario <- function(scenario, B = 1000, n_cores = 2, seed = 1234) {
  # scenario: list with elements n, true_hr, dist, median_surv, surv48, P_drop, beta, gamma, delta_baseline
  set.seed(seed)
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  # Export necessary functions and data to cluster
  clusterExport(cl, varlist = ls(envir = environment()), envir = environment())
  # Wrapper to run one replicate
  replicate_fun <- function(i) {
    run_single_rep(
      n = scenario$n,
      true_hr = scenario$true_hr,
      dist = scenario$dist,
      median_surv = scenario$median_surv,
      surv48 = scenario$surv48,
      P_drop = scenario$P_drop,
      beta = scenario$beta,
      gamma = scenario$gamma,
      delta_baseline = scenario$delta_baseline,
      prognostic_spec = scenario$prognostic_spec
    )
  }
  # Run B replicates in parallel
  reps <- parLapply(cl, seq_len(B), replicate_fun)
  # Combine results
  df <- do.call(rbind, reps)
  # Aggregate summary (Monte Carlo averages)
  summary_df <- aggregate(cbind(hr, bias, coverage, ESS_treated,
                                std_diff_X1, std_diff_X2, std_diff_X3, std_diff_X4) ~ method,
                          data = df, FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                         se = sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))))
  # Flatten columns
  summary_flat <- do.call(data.frame, summary_df)
  return(list(per_rep = df, summary = summary_flat))
}

# ------------------------
# Example: setting up and running a small grid
# ------------------------

# Example fixed beta vector (user should calibrate this to achieve desired control event rates separately)
# Here we assume moderate prognostic effects
beta_example <- c(log(1.3), log(1.3), log(1.3))  # effects for X1, X2, X3

# Calibrate delta for desired overlaps (e.g., 0.15, 0.25, 0.33, 0.5)
mu0 <- c(0.5, 0.5, 0.5, 0.5)
target_overlaps <- c(0.15, 0.25, 0.33, 0.5)
delta_list <- sapply(target_overlaps, function(p) {
  calibrate_delta(p, mu0 = mu0, prognostic_idx = 1:3, Nlarge = 10000)
})
names(delta_list) <- paste0("overlap_", target_overlaps)

# Calibrate gamma for a target ATUT HR (example for HR=0.7, Weibull)
# Warning: this is expensive; consider caching results
gamma_calibrated <- calibrate_gamma(target_hr = 0.7, beta = beta_example, dist = "weibull",
                                    median_surv = 24, surv48 = 0.05, P_drop = 0.2,
                                    target_overlap_delta = 1, Ncalib = 2000)  # reduced for example speed

# Define a scenario
scenario1 <- list(
  n = 200,
  true_hr = 0.7,
  dist = "weibull",
  median_surv = 24,
  surv48 = 0.05,
  P_drop = 0.2,
  beta = beta_example,
  gamma = gamma_calibrated,
  delta_baseline = delta_list["overlap_0.25"],  # e.g., target 25% overlap
  prognostic_spec = "correct"  # can be "correct", "under", "over"
)

# Run small-scale simulation (e.g., 100 replicates, 2 cores)
res <- run_scenario(scenario1, B = 100, n_cores = 2, seed = 2025)

# Inspect summary
print(res$summary)

# ------------------------
# Notes:
# - Calibration of beta to target control event rate is not fully automated here; one would apply a similar root-finding loop
#   varying a scaling on beta (or baseline hazard) to hit the desired marginal event rate (accounting for dropout).
# - MEW variants are placeholders; implementing maximum-ESS weights under exact moment constraints requires quadratic programming
#   (minimizing sum(w^2) subject to linear constraints) and is left for an extended implementation.
# - The ATUT HR is estimated by reweighting treated to control prognostic covariate distribution and fitting a weighted Cox model.
# - Robustness scenarios ("under", "over") control inclusion/omission of prognostic and nuisance covariates in the weighting step.
# - The structure allows extension: grid expansion over distributions, HR targets, dropout rates, sample sizes, overlap levels, etc.
# - Aggregation of Monte Carlo standard errors is done in summary; further bootstrap or variance stabilization can be added as needed.
