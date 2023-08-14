# Functions for adjusting step for anchored cases

# functions to be exported ---------------------------------------

#' Bucher method for adjusted treatment effect
#'
#' Given two estimated treatment effects of A vs. C and B vs. C
#' (i.e. two point estimates and corresponding standard errors),
#' derive the adjusted treatment of A vs. B using Bucher method,
#' with two-sided confidence limits and Z-test p-value
#'
#' @param trt a list with two named scalars for the study with interested experimental arm,
#' one named `'est'` for the point estimate, and the other named `'se'` for the standard error
#' @param com same as \code{trt}, but for the study with interested control arm
#' @param conf_lv a numerical scalar, prescribe confidence level to derive two-sided
#' confidence interval for the adjusted treatment effect
#'
#' @return a list with 5 elements,
#' \describe{
#'   \item{est}{a scalar, point estimate of the adjusted treatment effect}
#'   \item{se}{a scalar, standard error of the adjusted treatment effect (i.e. \code{est} in return)}
#'   \item{ci_l}{a scalar, lower confidence limit of a two-sided CI with prescribed nominal level by \code{conf_lv}}
#'   \item{ci_u}{a scalar, upper confidence limit of a two-sided CI with prescribed nominal level by \code{conf_lv}}
#'   \item{pval}{p-value of Z-test, with null hypothesis that \code{est} is zero}
#' }
#' @examples
#' @export

bucher <- function(trt, com, conf_lv = 0.95) {
  est <- trt$est - com$est
  se <- sqrt(trt$se^2 + com$se^2)
  ci_l <- est - stats::qnorm(0.5 + conf_lv / 2) * se
  ci_u <- est + stats::qnorm(0.5 + conf_lv / 2) * se
  if (est > 0) {
    pval <- 2 * (1 - stats::pnorm(est, 0, se))
  } else {
    pval <- 2 * stats::pnorm(est, 0, se)
  }

  list(
    est = est,
    se = se,
    ci_l = ci_l,
    ci_u = ci_u,
    pval = pval
  )
}

# functions NOT to be exported ---------------------------------------

#' Report-friendly output format for result from Bucher's method
#'
#' @param output output from \code{\link{bucher}} function
#' @param ci_digits an integer, number of decimal places for point estimate and derived confidence limits
#' @param pval_digits an integer, number of decimal places to display Z-test p-value
#'
#' @return a character vector of two elements, first element in format of 'est (ci_l; ci_u)'`,
#'  second element is the Z-test p-value, rounded according to \code{pval_digits}

print_bucher <- function(output, ci_digits = 2, pval_digits = 3) {
  res <- paste0(
    format(round(output$est, ci_digits), nsmall = ci_digits), " (",
    format(round(output$ci_l, ci_digits), nsmall = ci_digits), ";",
    format(round(output$ci_u, ci_digits), nsmall = ci_digits), ")"
  )

  disp_pval <- round(output$pval, pval_digits)
  disp_pval <- ifelse(disp_pval == 0, paste0("<", 1 / (10^pval_digits)), format(output$pval, nsmall = pval_digits))

  c(result = res, pvalue = disp_pval)
}
