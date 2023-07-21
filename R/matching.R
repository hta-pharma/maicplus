# Functions for matching step: estimation of individual weights

# functions to be exported ---------------------------------------

#' Derive individual weights in the matching step of MAIC
#'
#' Assuming data is properly processed, this function takes individual patient data (IPD) with centered covariates
#' (effect modifiers and/or prognostic variables) as input, and generates weights for each individual in IPD trial that
#' matches the chosen statistics of those covariates in Aggregated Data (AgD) trial.
#'
#' @param data a numeric matrix, centered covariates of IPD, no missing value in any cell is allowed
#' @param centered_colnames a character or numeric vector (column indicators) of centered covariates
#' @param start_val a scalar, the starting value for all coefficients of the propensity score regression
#' @param method a string, name of the optimization algorithm (see 'method' argument of \code{base::optim()}).
#' The default is `"BFGS"`, other options are `"Nelder-Mead"`, `"CG"`, `"L-BFGS-B"`, `"SANN"`, and `"Brent"`
#' @param ... all other arguments from \code{base::optim()}
#'
#' @return a list with the following 4 elements,
#' \describe{
#'   \item{data}{a data.frame, includes the input \code{data} with appended column 'weights' and 'scaled_weights'.
#'   Scaled weights has a summation to be the number of rows in \code{data} that has no missing value in any of the
#'   effect modifiers}
#'   \item{centered_colnames}{column names of centered effect modifiers in \code{data}}
#'   \item{nr_missing}{number of rows in \code{data} that has at least 1 missing value in specified centered effect
#'   modifiers}
#'   \item{ess}{effective sample size, square of sum divided by sum of squares}
#'   \item{opt}{R object returned by \code{base::optim()}, for assess convergence and other details}
#' }
#' @export
#'
#' @examples
estimate_weights <- function(data, centered_colnames = NULL, start_val = 0, method = "BFGS", ...) {
  # pre check
  ch1 <- is.data.frame(data)
  if (!ch1) stop("'data' is not a data.frame")

  ch2 <- (!is.null(centered_colnames))
  if (ch2 && is.numeric(centered_colnames)) {
    ch2b <- any(centered_colnames < 1 | centered_colnames > ncol(data))
    if (ch2b) stop("specified centered_colnames are out of bound")
  } else if (ch2 && is.character(centered_colnames)) {
    ch2b <- !all(centered_colnames %in% names(data))
    if (ch2b) stop("1 or more specified centered_colnames are not found in 'data'")
  } else {
    stop("'centered_colnames' should be either a numeric or character vector")
  }

  ch3 <- sapply(seq_along(centered_colnames), function(ii) {
    !is.numeric(data[, centered_colnames[ii]])
  })
  if (any(ch3)) {
    stop(paste0("following columns of 'data' are not numeric for the calculation:", paste(which(ch3), collapse = ",")))
  }

  # prepare data for optimization
  if (is.null(centered_colnames)) centered_colnames <- seq_len(ncol(data))
  EM <- data[, centered_colnames, drop = FALSE]
  ind <- apply(EM, 1, function(xx) any(is.na(xx)))
  nr_missing <- sum(ind)
  EM <- as.matrix(EM[!ind, , drop = FALSE])

  # objective and gradient functions
  objfn <- function(alpha, X) {
    sum(exp(X %*% alpha))
  }
  gradfn <- function(alpha, X) {
    colSums(sweep(X, 1, exp(X %*% alpha), "*"))
  }

  # estimate weights
  opt1 <- optim(
    par = rep(start_val, ncol(EM)),
    fn = objfn, gr = gradfn,
    X = EM,
    method = method,
    control = list(maxit = 300, trace = 2), ...
  )
  alpha <- opt1$par
  wt <- exp(EM %*% alpha)
  wt_rs <- (wt / sum(wt)) * nrow(EM)

  # append weights to data
  data$weights <- NA
  data$weights[!ind] <- wt

  data$scaled_weights <- NA
  data$scaled_weights[!ind] <- wt_rs

  if (is.numeric(centered_colnames)) centered_colnames <- names(data)[centered_colnames]

  # Output
  list(
    data = data,
    centered_colnames = centered_colnames,
    nr_missing = nr_missing,
    ess = sum(wt)^2 / sum(wt^2),
    opt = opt1
  )
}


#' Plot MAIC weights in a histogram with key statistics in legend
#'
#' Generates a plot given the individuals weights with key summary in top right legend that includes
#' median weight, effective sample size (ESS), and reduction percentage (what percent ESS is reduced from the
#' original sample size).
#' There are two options of weights provided in \code{\link{estimate_weights}}: unscaled or scaled.
#' Scaled weights are relative to the original unit weights of each individual.
#' In other words, a scaled weight greater than 1 means that an individual carries more weight in the
#' re-weighted population than the original data and a scaled weight less than 1 means that an individual carries
#' less weight in the re-weighted population than the original data.
#'
#' @param wt a numeric vector of individual MAIC weights (derived using \code{\link{estimate_weights}})
#' @param main_title a character string, main title of the plot
#'
#' @return a plot of unscaled or scaled weights
#' @importFrom graphics hist
#' @export

plot_weights <- function(wt, main_title = "Unscaled Individual Weights") {
  # calculate sample size and exclude NA from wt
  nr_na <- sum(is.na(wt))
  n <- length(wt) - nr_na
  wt <- na.omit(wt)

  # calculate effective sample size (ESS) and reduction of original sample
  ess <- (sum(wt)^2) / sum(wt^2)
  ess_reduct <- round((1 - (ess / n)) * 100, 2)

  # prepare legend
  plot_legend <- c(
    paste0("Median = ", round(median(wt), 4)),
    paste0("ESS = ", round(ess, 2)),
    paste0("Reduction% = ", ess_reduct)
  )
  plot_lty <- c(2, NA, NA)

  if (nr_na > 0) {
    plot_legend <- c(plot_legend, paste0("#Missing Weights = ", nr_na))
    plot_lty <- c(plot_lty, NA)
  }

  # plot
  par(mfrow = c(1, 1), family = "HersheySerif", mgp = c(2.3, 0.5, 0), cex.axis = 0.9, cex.lab = 0.95, bty = "n")
  hist(wt, border = "white", col = "#6ECEB2", main = main_title, breaks = 20, yaxt = "n", xlab = "")
  axis(2, las = 1)
  abline(v = median(wt), lty = 2, col = "#688CE8", lwd = 2)
  legend("topright", bty = "n", lty = plot_lty, cex = 0.8, legend = plot_legend)
}

#' Check to see if weights are optimized correctly
#'
#' This function checks to see if the optimization is done properly by checking the covariate averages
#' before and after adjustment.
#'
#' @param optimized object returned after calculating weights using \code{\link{estimate_weights}}
#' @param match_cov covariates that should be checked to see if the IPD weighted average matches the aggregate data
#' average. This could be same set of variables that were used to match or it can include variables that were not
#' included to match (i.e. stratification variables)
#' @param digits number of digits for rounding summary table
#'
#' @return data.frame of weighted and unweighted covariate averages of the IPD
#' @export



check_weights <- function(optimized, match_cov, digits = 2) {
  if (missing(match_cov)) stop("match_cov is missing. Covariates to check must be defined.")
  ipd_with_weights <- optimized$data

  arm_names <- c("Unweighted IPD", "Weighted IPD")
  ess <- c(nrow(ipd_with_weights), optimized$ess)

  unweighted_cov <- sapply(ipd_with_weights[, match_cov, drop = FALSE], mean)

  weighted_cov <- sapply(
    ipd_with_weights[, match_cov, drop = FALSE],
    weighted.mean,
    w = ipd_with_weights$weights
  )

  cov_combined <- rbind(unweighted_cov, weighted_cov)

  baseline_summary <- cbind(ess, cov_combined)
  rownames(baseline_summary) <- arm_names

  round(baseline_summary, digits = digits)
}