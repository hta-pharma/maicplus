#' Bucher method for adjusted treatment effect
#'
#' Given two estimated treatment effects of A vs. C and B vs. C (i.e. two point estimates and corresponding standard errors),
#' derive the adjusted treatment of A vs. B using Bucher method, with two-sided confidence limits and Z-test p-value
#'
#' @param trt a list with two named scalars for the study with interested experimental arm, one named 'est' for the point estimate, and the other named 'se' for the standard error
#' @param com same as \code{trt} for the study with interested control arm
#' @param conf.lv a numerical scalar, prescribe confidence level to derive two-sided confidence interval for the adjusted treatment effect
#'
#' @return a list with 5 elements,
#' \describe{
#'   \item est - a scalar, point estimate of the adjusted treatment effect
#'   \item se - a scalar, standard error of the adjusted treatment effect (i.e. \code{est} in return)
#'   \item ci_l - a scalor, lower confidence limit of a two-sided CI with prescribed nominal level by \code{conf.lv}
#'   \item ci_u - a scalor, upper confidence limit of a two-sided CI with prescribed nominal level by \code{conf.lv}
#'   \item pval - p-value of Z-test, with null hypothesis that \code{est} is zero
#' }
#' @export
#'
#' @examples
bucher <- function(trt, com, conf.lv = 0.95) {
  est <- trt$est - com$est
  se <- sqrt(trt$se^2 + com$se^2)
  ci_l <- est - qnorm(0.5 + conf.lv / 2) * se
  ci_u <- est + qnorm(0.5 + conf.lv / 2) * se
  if (est > 0) {
    pval <- 2 * (1 - pnorm(est, 0, se))
  } else {
    pval <- 2 * pnorm(est, 0, se)
  }

  list(
    est = est,
    se = se,
    ci_l = ci_l,
    ci_u = ci_u,
    pval = pval
  )
}

#' Report-friendly output format for resutl from Bucher's method
#'
#' @param output output from \code{\link{bucher}} function
#' @param ci.digits an integer, number of decimal places for point estimate and derived confidence limits
#' @param pval.digits an integer, number of decimal places to display Z-test p-value
#'
#' @return a character vector of two elements, first element inf format of 'est (ci_l; ci_u)', second element is Z-test p-value, rounded according to \code{pval.digits}

print.bucher <- function(output, ci.digits = 2, pval.digits = 3) {
  res <- paste0(
    format(round(output$est, ci.digits), nsmall = ci.digits), " (",
    format(round(output$ci_l, ci.digits), nsmall = ci.digits), ";",
    format(round(output$ci_u, ci.digits), nsmall = ci.digits), ")"
  )

  disp.pval <- round(output$pval, pval.digits)
  disp.pval <- ifelse(disp.pval == 0, paste0("<", 1 / (10^pval.digits)), format(output$pval, nsmall = pval.digits))

  c(result = res, pvalue = disp.pval)
}
