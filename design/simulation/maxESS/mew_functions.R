

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