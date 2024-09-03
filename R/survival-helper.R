#' Helper function to retrieve median survival time from a `survival::survfit` object
#'
#' Extract and display median survival time with confidence interval
#'
#' @param km_fit returned object from \code{survival::survfit}
#' @param legend a character string, name used in 'type' column in returned data frame
#' @param time_scale a character string, 'years', 'months', 'weeks' or 'days', time unit of median survival time
#'
#' @examples
#' data(adtte_sat)
#' data(pseudo_ipd_sat)
#' library(survival)
#' combined_data <- rbind(adtte_sat[, c("TIME", "EVENT", "ARM")], pseudo_ipd_sat)
#' kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data, conf.type = "log-log")
#'
#' # Derive median survival time
#' medSurv <- medSurv_makeup(kmobj, legend = "before matching", time_scale = "day")
#' medSurv
#' @return a data frame with a index column 'type', median survival time and confidence interval
#' @export

medSurv_makeup <- function(km_fit, legend = "before matching", time_scale) {
  time_scale <- match.arg(time_scale, choices = c("years", "months", "weeks", "days"))

  # km_fit is the returned object from survival::survfit
  km_fit <- summary(km_fit)$table
  km_fit[, 5:ncol(km_fit)] <- get_time_as(km_fit[, 5:ncol(km_fit)], time_scale)

  toyadd <- data.frame(
    treatment = gsub("ARM=", "", rownames(km_fit)),
    type = rep(legend, 2)
  )

  km_fit <- cbind(toyadd, km_fit)
  rownames(km_fit) <- NULL

  km_fit
}


#' Helper function to select set of variables used for Kaplan-Meier plot
#'
#' @param km_fit returned object from \code{survival::survfit}
#' @param single_trt_name name of treatment if no strata are specified in `km_fit`
#'
#' @examples
#' library(survival)
#' data(adtte_sat)
#' data(pseudo_ipd_sat)
#' combined_data <- rbind(adtte_sat[, c("TIME", "EVENT", "ARM")], pseudo_ipd_sat)
#' kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data, conf.type = "log-log")
#' survfit_makeup(kmobj)
#' @return a list of data frames of variables from [survival::survfit()]. Data frame is divided by treatment.
#' @export

survfit_makeup <- function(km_fit, single_trt_name = "treatment") {
  # in case km_fit is only for single arm
  if ("strata" %in% names(km_fit)) {
    use_trt <- mapply(rep, 1:2, each = km_fit$strata)
    if (is.list(use_trt)) use_trt <- unlist(use_trt)
    if (is.matrix(use_trt)) use_trt <- as.vector(use_trt)
    is_single <- FALSE
  } else {
    use_trt <- rep(single_trt_name, length(km_fit$time))
    is_single <- TRUE
  }

  kmdat <- data.frame(
    time = km_fit$time,
    treatment = use_trt,
    n.risk = km_fit$n.risk,
    n.event = km_fit$n.event,
    censor = km_fit$n.censor,
    surv = km_fit$surv,
    lower = km_fit$lower,
    upper = km_fit$upper,
    cumhaz = km_fit$cumhaz
  )
  if (!is_single) kmdat$treatment <- sapply(strsplit(names(km_fit$strata), "="), "[[", 2)[kmdat$treatment]
  split(kmdat, f = kmdat$treatment)
}
