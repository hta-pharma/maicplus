
#' Calculate standard error from the reported confidence interval.
#'
#' Often times in analysis of clinical trials, standard error is not reported and only the confidence interval is reported.
#' For instance, a given study can report a hazard ratio of 0.70 and 98.22% confidence interval of 0.55 to 0.90.
#' This function calculates standard error of the estimate given the reported confidence interval.
#'  
#' @param RR Reported estimate of relative risk (i.e. hazard ratio or odds ratio). Note that this estimate is not logged.
#' @param CI_lower Reported lower percentile value of the relative risk
#' @param CI_upper Reported upper percentile value of the relative risk 
#' @param CI_perc Reported percentage of confidence interval reported
#' @return Standard error of the relative risk and log relative risk
#' @examples
#' find_SE_fromCI(RR = 0.70, CI_lower = 0.55, CI_upper = 0.90, CI_perc = 0.95)
#' @export

find_SE_fromCI <- function(RR = NULL, CI_lower = NULL, CI_upper = NULL, 
                           CI_perc = 0.95){
  
  logRR <- log(RR)
  alpha <- 1 - CI_perc
  logRR_SE <- (log(CI_upper) - log(CI_lower)) / (2 * qnorm(1 - alpha/2))
  RR_SE <- logRR_SE * RR
  
  return(list(logRR_SE = logRR_SE, RR_SE = RR_SE))
}

#' Retrieve relative risk and log relative risk with associated standard error and confidence interval
#' 
#' Convenient function to retrieve both the relative risk (RR) and log relative risk (logRR) either from
#' fitted model or separate calculation.
#' Relative risk would be hazard ratio for the cox regression and odds ratio for the logistic regression.
#' Delta method is used to calculate the standard error of the relative risk. 
#' There are two ways of finding CI of the relative risk: 
#' use se(RR) obtained via the delta method and calculate the end points or 
#' form a confidence interval for logRR and then exponentiate the endpoints. 
#' We will use the second method since logRR converges to a normal distribution 
#' more quickly than RR. This is also a default coxph method in [survival] R package.
#'  
#' @param fit weighted or unweighted regression fit (i.e. cox regression or logistic regression) where the logRR can be extracted from
#' @param logRR Instead of providing the fitted model (via parameter [fit]), one can specify logRR manually through this parameter
#' @param logRR_SE If logRR is specified, standard error has to be provided.
#' @param conf_lv Defines the percentage of confidence interval 
#' @return a list with 6 elements,
#' \describe{
#'   \item{RR}{an estimate of relative risk}
#'   \item{RR_SE}{standard error of the relative risk}
#'   \item{RR_CI}{Two-sided CI with prescribed nominal level by [conf_lv]}
#'   \item{pval}{p-value of Z-test, with null hypothesis that RR is zero}
#'   \item{logRR}{an estimate of log of relative risk}
#'   \item{logRR_SE}{standard error of the log of relative risk}
#' }
#' @example inst/examples/find_RR_ex.R
#' @export

find_RR <- function(fit = NULL, logRR = NULL, logRR_SE = NULL, conf_lv = 0.95){
  
  if(!is.null(fit)){
    if(inherits(fit, "coxph")){
      if(is.null(fit$weights)){
        logRR <- as.vector(coef(fit))
        logRR_SE <- summary(fit)$coefficients[,"se(coef)"]
      } else{
        logRR <- as.vector(coef(fit))
        logRR_SE <- summary(fit)$coefficients[,"robust se"]
      }
    } else if(inherits(fit, "glm")){
      if(length(unique(fit$prior.weights)) == 1){
        logRR <- as.vector(coef(fit)[-1])
        logRR_SE <- summary(fit)$coefficients[-1, "Std. Error"]  
      } else{
        logRR <- as.vector(coef(fit)[-1])
        V.sw <- sandwich::vcovHC(fit, type = "HC3")
        logRR_SE <- sqrt(V.sw[2, 2])
      }
    }
  } else {
    if(is.null(logRR) || is.null(logRR_SE)){
      stop("If fitted model has not be specified, logRR and logRR_SE have to be specified")
    }
  }
  
  RR <- as.vector(exp(logRR))
  RR_SE <- as.vector(logRR_SE * RR)
  
  RR_CI <- exp(c(logRR - qnorm(1 - (1-conf_lv)/2) * logRR_SE, 
                 logRR + qnorm(1 - (1-conf_lv)/2) * logRR_SE))
  names(RR_CI) <- c("lower", "upper")
  
  if (RR > 0) {
    pval <- 2 * (1 - stats::pnorm(RR, 0, RR_SE))
  } else {
    pval <- 2 * stats::pnorm(RR, 0, RR_SE)
  }
  
  return(list(RR = RR, RR_SE = RR_SE, RR_CI = RR_CI, logRR = logRR, pval = pval, logRR_SE = logRR_SE))
}


#' Calculate relative risk using bucher's method for anchored comparison
#'
#' Convenient function to calculate relative risk using bucher's method.
#' Function calculates treatment effect of C vs B via subtracting treatment effect of C vs A from treatment effect of B vs A.
#' Function also calculates standard error and confidence interval for the relative risk (C vs B).
#' This function is for an anchored comparison as it assumes a common comparator "A".
#' Delta method is used to calculate the standard error. 
#' There are two ways of finding CI of the relative risk: 
#' use se(RR) obtained via the delta method and calculate the end points or 
#' form a confidence interval for logRR and then exponentiate the endpoints. 
#' We will use the second method since logRR converges to a normal distribution 
#' more quickly than RR. This is also a default coxph method in [survival] R package.
#'  
#' @param trt A list of two scalars for the study with interested experimental arm.
#' One named `'logRR'` for the log relative risk and the other named `'logRR_SE'` for the standard error of log relative risk
#' One can also use the object returned from [find_RR] for convenience.
#' If only `'RR'` and `'RR_SE'` are available for the internal study, function
#' first calculates `'logRR'` and `'logRR_SE'` from the provided relative risk.
#' @param com Same list as [trt], but for the comparator study
#' If only `'RR'` and `'RR_SE'` are available for the comparator study, function
#' first calculates `'logRR'` and `'logRR_SE'` from the provided relative risk.
#' @param conf_lv Defines the percentage of confidence interval
#' @return a list with 6 elements,
#' \describe{
#'   \item{RR}{an estimate of relative risk}
#'   \item{RR_SE}{standard error of the relative risk}
#'   \item{RR_CI}{Two-sided CI with prescribed nominal level by [conf_lv]}
#'   \item{pval}{p-value of Z-test, with null hypothesis that RR is zero}
#'   \item{logRR}{an estimate of log of relative risk}
#'   \item{logRR_SE}{standard error of the log of relative risk}
#' }
#' @example inst/examples/bucher_ex.R
#' @export

bucher <- function(trt = NULL, com = NULL, conf_lv = 0.95){
  
  if(is.null(trt$logRR) & !is.null(trt$RR) & !is.null(trt$RR_SE)){
    trt$logRR <- log(trt$RR)
    trt$logRR_SE <- trt$RR_SE/trt$RR
  }
  
  if(is.null(com$logRR) & !is.null(com$RR) & !is.null(com$RR_SE)){
    com$logRR <- log(com$RR)
    com$logRR_SE <- com$RR_SE/com$RR
  }
  
  logRR <- trt$logRR - com$logRR
  logRR_SE <- sqrt(trt$logRR_SE^2 + com$logRR_SE^2)
  
  RR <- exp(logRR)
  RR_SE <- logRR_SE * RR
  RR_CI <- exp(c(logRR - qnorm(1 - (1-conf_lv)/2) * logRR_SE, 
                 logRR + qnorm(1 - (1-conf_lv)/2) * logRR_SE))
  names(RR_CI) <- c("lower", "upper")
  
  if (RR > 0) {
    pval <- 2 * (1 - stats::pnorm(RR, 0, RR_SE))
  } else {
    pval <- 2 * stats::pnorm(RR, 0, RR_SE)
  }
  
  outdata <- list(RR = RR, RR_SE = RR_SE, RR_CI = RR_CI, pval = pval, logRR = logRR, logRR_SE = logRR_SE)
  
  class(outdata) <- c("maicplus_bucher", "list")
  outdata
  
}


#' Print method for bucher object
#'
#' @param x object from [bucher]
#' @param ci_digits an integer, number of decimal places for point estimate and derived confidence limits
#' @param pval_digits an integer, number of decimal places to display Z-test p-value
#' @describeIn bucher Print method for bucher objects
#' @export

print.maicplus_bucher <- function(x, ci_digits = 2, pval_digits = 3) {
  res <- paste0(
    format(round(x$RR, ci_digits), nsmall = ci_digits), " [",
    format(round(x$RR_CI[1], ci_digits), nsmall = ci_digits), ";",
    format(round(x$RR_CI[2], ci_digits), nsmall = ci_digits), "]"
  )
  
  disp_pval <- round(x$pval, pval_digits)
  disp_pval <- ifelse(disp_pval == 0, paste0("<", 1 / (10^pval_digits)), format(disp_pval, nsmall = pval_digits))
  
  output <- c(res, disp_pval)
  names(output) <- c("result", "pvalue")
  return(output)
}
