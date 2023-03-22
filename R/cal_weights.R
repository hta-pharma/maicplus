#' Derive individual weights in the matching step of MAIC
#'
#' @param EM a numeric matrix, centered effective modifiers of IPD, no missing value in any cell is allowed
#' @param startVal a scalar, the starting value for all coefficients of the propensity score regression
#' @param ... all other arguments from \code{base::optim()}
#'
#' @return a list with 4 elements,
#' \describe{
#'   \item wt - a numeric vector of unscaled individual weights.
#'   \item wt.rs - a numerical vector of rescaled individual weights, with summation equaling to sample size (# rows of input \code{EM})
#'   \item ess - effective sample size, square of sum divided by sum of squares
#'   \item opt - R object returned by \code{base::optim()}, for assess convergence and other details
#' }
#' @export
#'
#' @examples
cal_weights <- function(EM, startVal = 0, method = "BFGS", ...) {
  ch1 <- apply(EM, 1, function(xxx) any(is.na(xxx)))
  if (any(ch1)) stop(paste0("Following rows has missing value: ", paste(which(ch1), collapse = ",")))

  objfn <- function(alpha, X) {
    sum(exp(X %*% alpha))
  }
  gradfn <- function(alpha, X) {
    colSums(sweep(X, 1, exp(X %*% alpha), "*"))
  }

  # Estimate weights
  opt1 <- optim(
    par = rep(startVal, ncol(EM)),
    fn = objfn, gr = gradfn,
    X = EM,
    method = method,
    control = list(maxit = 300, trace = 2), ...
  )
  alpha <- opt1$par
  wt <- exp(EM %*% alpha)
  wt.rs <- (wt / sum(wt)) * nrow(EM)

  # Output
  list(
    wt = wt,
    wt.rs = wt.rs,
    ess = sum(wt)^2 / sum(wt^2),
    opt = opt1
  )
}
