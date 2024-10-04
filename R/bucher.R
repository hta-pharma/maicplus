#' Bucher method for combining treatment effects
#'
#' Given two treatment effects of A vs. C and B vs. C
#' derive the treatment effects of A vs. B using the Bucher method.
#' Two-sided confidence interval and Z-test p-value are also calculated.
#' Treatment effects and standard errors should be in log scale
#' for hazard ratio, odds ratio, and risk ratio.
#' Treatment effects and standard errors should be in natural scale
#' for risk difference and mean difference.
#'
#' @param trt a list of two scalars for the study with the
#' experimental arm. `'est'` is the point estimate and `'se'`
#' is the standard error of the treatment effect.
#' For time-to-event data, `'est'` and `'se'` should be point estimate and
#' standard error of the log hazard ratio.
#' For binary data, `'est'` and `'se'` should be point estimate and
#' standard error of the log odds ratio, log risk ratio, or risk
#' difference.
#' For continuous data,  `'est'` and `'se'` should be point estimate and
#' standard error of the mean difference.
#' @param com same as \code{trt}, but for the study with the
#' control arm
#' @param conf_lv a numerical scalar, prescribe confidence level to derive
#' two-sided confidence interval for the treatment effect
#'
#' @return a list with 5 elements,
#' \describe{
#'   \item{est}{a scalar, point estimate of the treatment effect}
#'   \item{se}{a scalar, standard error of the treatment effect}
#'   \item{ci_l}{a scalar, lower confidence limit of a two-sided CI
#'   with prescribed nominal level by \code{conf_lv}}
#'   \item{ci_u}{a scalar, upper confidence limit of a two-sided CI
#'   with prescribed nominal level by \code{conf_lv}}
#'   \item{pval}{p-value of Z-test, with null hypothesis that
#'   \code{est} is zero}
#' }
#' @export
#' @examples
#' trt <- list(est = log(1.1), se = 0.2)
#' com <- list(est = log(1.3), se = 0.18)
#' result <- bucher(trt, com, conf_lv = 0.9)
#' print(result, ci_digits = 3, pval_digits = 3)
bucher <- function(trt, com, conf_lv = 0.95) {
  if (!isTRUE(is.finite(trt$est))) stop("trt$est is not valid: ", trt$est)
  if (!isTRUE(is.finite(trt$se))) stop("trt$se is not valid: ", trt$se)
  if (!isTRUE(is.finite(com$est))) stop("com$est is not valid: ", com$est)
  if (!isTRUE(is.finite(com$se))) stop("com$se is not valid: ", com$se)
  if (conf_lv < 0 || 1 < conf_lv) stop("conf_lv must be in (0, 1): ", conf_lv)

  est <- trt$est - com$est
  se <- sqrt(trt$se^2 + com$se^2)
  ci_l <- est - stats::qnorm(0.5 + conf_lv / 2) * se
  ci_u <- est + stats::qnorm(0.5 + conf_lv / 2) * se
  if (est > 0) {
    pval <- 2 * (1 - stats::pnorm(est, 0, se))
  } else {
    pval <- 2 * stats::pnorm(est, 0, se)
  }

  outdata <- list(
    est = est,
    se = se,
    ci_l = ci_l,
    ci_u = ci_u,
    pval = pval
  )

  class(outdata) <- c("maicplus_bucher", "list")
  outdata
}


#' Calculate standard error from the reported confidence interval.
#'
#' Comparator studies often only report confidence interval of the
#' treatment effects. This function calculates standard error of the
#' treatment effect given the reported confidence interval.
#' For relative treatment effect (i.e. hazard ratio, odds ratio, and
#' risk ratio), the function would log the confidence interval.
#' For risk difference and mean difference,
#' we do not log the confidence interval.
#' The option to log the confidence interval is controlled
#' by `'log'` parameter.
#'
#' @param CI_lower Reported lower percentile value of the
#' treatment effect
#' @param CI_upper Reported upper percentile value of the
#' treatment effect
#' @param CI_perc Percentage of confidence interval reported
#' @param log Whether the confidence interval should be logged.
#' For relative treatment effect, log should be applied because
#' estimated log treatment effect is approximately normally distributed.
#' @return Standard error of log relative treatment effect if `'log'`
#' is true and standard error of the treatment effect if `'log'`
#' is false
#' @examples
#' find_SE_from_CI(CI_lower = 0.55, CI_upper = 0.90, CI_perc = 0.95)
#' @export

find_SE_from_CI <- function(CI_lower = NULL, CI_upper = NULL,
                            CI_perc = 0.95, log = TRUE) {
  if (CI_perc > 1 || CI_perc < 0) {
    stop("CI_perc has to be between 0 and 1")
  }

  if (is.null(CI_lower) || is.null(CI_upper)) {
    stop("Both CI_lower and CI_upper need to be specified")
  }

  if (!is.numeric(CI_lower) || !is.numeric(CI_upper)) {
    stop("Both CI_lower and CI_upper need to be specified")
  }

  alpha <- 1 - CI_perc
  se <- ifelse(log,
    (log(CI_upper) - log(CI_lower)) / (2 * qnorm(1 - alpha / 2)),
    (CI_upper - CI_lower) / (2 * qnorm(1 - alpha / 2))
  )
  return(se)
}

#' Print method for `maicplus_bucher` object
#'
#' @param x `maicplus_bucher` object
#' @param ci_digits an integer, number of decimal places for point
#' estimate and derived confidence limits
#' @param pval_digits an integer, number of decimal places to display
#' Z-test p-value
#' @param exponentiate whether the treatment effect and confidence
#' interval should be exponentiated. This applies to relative
#' treatment effects. Default is set to false.
#' @param ... not used
#' @describeIn bucher Print method for `maicplus_bucher` objects
#' @export

print.maicplus_bucher <- function(x, ci_digits = 2, pval_digits = 3,
                                  exponentiate = FALSE, ...) {
  output <- reformat(x, ci_digits, pval_digits,
    show_pval = TRUE, exponentiate
  )
  print(output)
}

#' Reformat `maicplus_bucher` alike object
#'
#' @param x a list, structured like a `maicplus_bucher` object
#' @param ci_digits an integer, number of decimal places for point
#' estimate and derived confidence limits
#' @param pval_digits an integer, number of decimal places to display
#' Z-test p-value
#' @param show_pval a logical value, default is TRUE. If FALSE, p-value will not
#' be output as the second element of the character vector
#' @param exponentiate whether the treatment effect and confidence
#' interval should be exponentiated. This applies to relative
#' treatment effects. Default is set to false.
#' @keywords internal

reformat <- function(x, ci_digits = 2, pval_digits = 3,
                     show_pval = TRUE, exponentiate = FALSE) {
  transform_this <- function(x) {
    ifelse(exponentiate, exp(x), x)
  }

  a <- format(round(transform_this(x$est), ci_digits),
    nsmall = ci_digits
  )
  b <- format(round(transform_this(x$ci_l), ci_digits),
    nsmall = ci_digits
  )
  c <- format(round(transform_this(x$ci_u), ci_digits),
    nsmall = ci_digits
  )
  res <- paste0(a, "[", b, "; ", c, "]")

  disp_pval <- round(x$pval, pval_digits)
  disp_pval <-
    ifelse(disp_pval == 0,
      paste0("<", 1 / (10^pval_digits)),
      format(disp_pval, nsmall = pval_digits)
    )

  if (show_pval) {
    output <- c(res, disp_pval)
    names(output) <- c("result", "pvalue")
  } else {
    output <- res
    names(output) <- "result"
  }

  output
}
