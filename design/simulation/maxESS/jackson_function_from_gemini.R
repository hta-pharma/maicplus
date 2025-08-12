#' Calculate MAIC weights that maximize effective sample size (ESS)
#'
#' This function computes MAIC weights by solving a constrained optimization
#' problem to maximize the effective sample size (ESS) while achieving exact
#' balance on baseline covariates. This methodology is an alternative to
#' conventional MAIC and is based on the work of Zubizarreta (2015), as
#' discussed in Jackson et al. (2020).
#'
#' @param covariates A data frame or matrix of individual-level baseline
#'   covariates from the trial of interest. Column names must match the
#'   names in `target_means`.
#' @param target_means A named numeric vector or list of the mean values of the
#'   same covariates in the target comparator population.
#'
#' @return A list containing:
#'   - `weights`: A numeric vector of the calculated MAIC weights.
#'   - `ess`: The effective sample size of the weighted population.
#'
#' @references
#' Jackson D, Rhodes K, Ouwens M. Alternative weighting schemes when
#' performing matching‐adjusted indirect comparisons. Res Syn Meth.
#' 2021;12:333-346. https://doi.org/10.1002/jrsm.1466
#'
#' @examples
#' # Simulate individual-level patient data from Trial A
#' set.seed(42)
#' n_a <- 100
#' trial_data <- data.frame(
#'   age = rnorm(n_a, mean = 65, sd = 8),
#'   sex = rbinom(n_a, 1, 0.6)
#' )
#'
#' # Define target covariate means from Comparator Trial B
#' target_means <- c(age = 68, sex = 0.5)
#'
#' # Ensure limSolve is installed before running
#' if (requireNamespace("limSolve", quietly = TRUE)) {
#'   # Calculate MAIC weights
#'   maic_result <- calculate_maic_jackson_weights(
#'     covariates = trial_data,
#'     target_means = target_means
#'   )
#'
#'   # Print the results
#'   print(maic_result)
#'
#'   # Check for balance: weighted means should be close to target means
#'   weighted_mean_age <- weighted.mean(trial_data$age, maic_result$weights)
#'   weighted_mean_sex <- weighted.mean(trial_data$sex, maic_result$weights)
#'
#'   cat("Target mean age:", target_means["age"], "| Weighted mean age:", round(weighted_mean_age, 3), "\n")
#'   cat("Target mean sex:", target_means["sex"], "| Weighted mean sex:", round(weighted_mean_sex, 3), "\n")
#' }
calculate_maic_jackson_weights <- function(covariates, target_means) {
  # Check for the required package
  if (!requireNamespace("limSolve", quietly = TRUE)) {
    stop("Package 'limSolve' needed for this function to work. Please install it.", call. = FALSE)
  }

  # Step 1: Center the covariates by subtracting the target population means.
  centered_covariates <- as.matrix(
    sweep(covariates[, names(target_means), drop = FALSE], 2, target_means, FUN = "-")
  )
  n <- nrow(centered_covariates)

  # Step 2: Set up the optimization problem for `limSolve::lsei`.
  # The objective is to minimize the sum of squared weights, min(||w*||^2).
  # We set the `A` matrix to the identity matrix and `B` vector to zeros.
  A <- diag(n)
  B <- rep(0, n)

  # The constraints are:
  # 1. Sum of weights equals 1
  # 2. Weighted mean of centered covariates equals zero
  # This corresponds to the `E` matrix and `F` vector arguments.
  E <- t(cbind(1, centered_covariates)) # The matrix from the paper (Z^T)
  F <- c(1, rep(0, ncol(centered_covariates)))

  # The inequality constraints require weights to be non-negative.
  # This corresponds to the `G` matrix and `H` vector.
  G <- diag(n)
  H <- rep(0, n)

  # Step 3: Solve the constrained least squares problem.
  # The `lsei` function finds the vector `w` that minimizes ||Aw - B||^2
  # subject to Ew = F and Gw >= H.
  optimization_result <- limSolve::lsei(
    A = A,
    B = B,
    E = E,
    F = F,
    G = G,
    H = H
  )

  # Step 4: Extract the weights from the result.
  weights <- optimization_result$X

  # Step 5: Calculate the effective sample size (ESS)
  # The paper's formula is (sum of weights)^2 / (sum of weights^2).
  # Since we constrained the weights to sum to 1, this simplifies to 1 / sum(weights^2).
  ess <- 1 / sum(weights^2)

  return(list(weights = weights, ess = ess))
}
