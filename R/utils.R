# Restore the RNG back to a previous state using the global .Random.seed
set_random_seed <- function(old_seed) {
  if (is.null(old_seed)) {
    rm(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    assign(".Random.seed", value = old_seed, envir = globalenv(), inherits = FALSE)
  }
}

construct_boot_data <- function(weighted_data, i = 1) {
  if (is.null(weighted_data$boot)) stop("Must contain bootstrap results from estimate_weights()")
  i <- as.integer(i)
  R <- dim(weighted_data$boot)[3]
  if (i < 1 || i > R) stop("i must be integer between 1 and ", R)

  boot_data <- weighted_data$boot[, , i]
  weighted_data$data <- weighted_data$data[boot_data[, 1], ]
  weighted_data$data$weights <- boot_data[, 2]
  weighted_data$data$scaled_weights <- boot_data[, 2] / sum(boot_data[, 2])
  weighted_data
}


transform_ratio <- function(object) {
  result <- object
  result$est <- exp(object$est)
  # log normal parameterization for SE
  result$se <- sqrt((exp(object$se^2) - 1) * exp(2 * object$est + object$se^2))
  result$ci_l <- exp(object$ci_l)
  result$ci_u <- exp(object$ci_u)
  result
}

transform_absolute <- function(object) {
  result <- object
  result$est <- object$est * 100
  result$se <- object$se * 100
  result$ci_l <- object$ci_l * 100
  result$ci_u <- object$ci_u * 100
  result
}
