# Functions for adjusting step for anchored cases

# functions to be exported ---------------------------------------

#' Bucher method for adjusted treatment effect
#'
#' Given two estimated treatment effects of A vs. C and B vs. C
#' (i.e. two point estimates and corresponding standard errors),
#' derive the adjusted treatment of A vs. B using Bucher method,
#' with two-sided confidence limits and Z-test p-value.
#' The relative treatment effects should be in log scale in the case of
#' hazard ratio, odds ratio, and risk ratio. The corresponding 
#' standard error would be standard error of the log relative treatment 
#' effect. 
#'
#' @param trt a list with two named scalars for the study with interested experimental arm,
#' one named `'est'` for the point estimate (log scale), and the other named `'se'` for the standard error
#' of the log treatment effect. For instance, for time-to-event data, `'est'` would be
#' log hazard ratio and `'se'` would be standard error of the log hazard ratio.
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
#' @export
#' @examples
#' trt <- list(est = log(1.1), se = 0.2)
#' com <- list(est = log(1.3), se = 0.18)
#' bucher(trt, com, conf_lv = 0.9)

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
#' risk ratio), the function would first log the confidence interval. 
#' For risk difference, we do not need to log the confidence interval.
#' The option to log the confidence interval is controlled 
#' by `'logged'` parameter.
#'  
#' @param CI_lower Reported lower percentile value of the relative risk
#' @param CI_upper Reported upper percentile value of the relative risk 
#' @param CI_perc Reported percentage of confidence interval reported
#' @return Standard error of log relative treatment effect if `'logged'` 
#' is true and standard error of the treatment effect if `'logged'`
#' is false
#' @examples
#' find_SE_fromCI(CI_lower = 0.55, CI_upper = 0.90, CI_perc = 0.95)
#' @export

find_SE_fromCI <- function(CI_lower = NULL, CI_upper = NULL, 
                           CI_perc = 0.95, logged = TRUE){
  
  alpha <- 1 - CI_perc
  
  if(logged){
    se <- (log(CI_upper) - log(CI_lower)) / (2 * qnorm(1 - alpha/2))
  } else{
    se <- (CI_upper - CI_lower) / (2 * qnorm(1 - alpha/2))
  }
  return(se)
}


# functions NOT to be exported ---------------------------------------

#' Print method for bucher object
#'
#' @param x object from [bucher]
#' @param ci_digits an integer, number of decimal places for point estimate and derived confidence limits
#' @param pval_digits an integer, number of decimal places to display Z-test p-value
#' @describeIn bucher Print method for bucher objects
#' @export

print.maicplus_bucher <- function(x, ci_digits = 2, pval_digits = 3) {
  res <- paste0(
    format(round(x$est, ci_digits), nsmall = ci_digits), "[",
    format(round(x$ci_l, ci_digits), nsmall = ci_digits), ";",
    format(round(x$ci_u, ci_digits), nsmall = ci_digits), "]"
  )

  disp_pval <- round(x$pval, pval_digits)
  disp_pval <- ifelse(disp_pval == 0, paste0("<", 1 / (10^pval_digits)), format(disp_pval, nsmall = pval_digits))

  output <- c(res, disp_pval)
  names(output) <- c("result", "pvalue")
  return(output)
}
