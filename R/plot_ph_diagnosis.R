#' Plot of Log Cumulative Hazard Rate
#'
#' a diagnosis plot for proportional hazard assumption, versus log-time (default) or time
#'
#' @param km_fit returned object from \code{survival::survfit}
#' @param time_scale a character string, 'year', 'month', 'week' or 'day', time unit of median survival time
#' @param log_time logical, TRUE (default) or FALSE
#' @param endpoint_name a character string, name of the endpoint
#' @param subtitle a character string, subtitle of the plot
#' @param exclude_censor logical, should censored data point be plotted
#' @examples
#' library(survival)
#' load(system.file("extdata", "combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
#' kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte, conf.type = "log-log")
#' log_cum_haz_plot(kmobj,
#'   time_scale = "month", log_time = TRUE,
#'   endpoint_name = "OS", subtitle = "(Before Matching)"
#' )
#' @return a plot of log cumulative hazard rate
#' @export

log_cum_haz_plot <- function(km_fit,
                             time_scale,
                             log_time = TRUE,
                             endpoint_name = "",
                             subtitle = "",
                             exclude_censor = TRUE) {
  time_unit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)
  if (!time_scale %in% names(time_unit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")

  clldat <- survfit_makeup(km_fit)

  if (exclude_censor) {
    clldat <- lapply(clldat, function(xxt) xxt[xxt$censor == 0, , drop = FALSE])
  }

  all.times <- do.call(rbind, clldat)$time / time_unit[[time_scale]]
  if (log_time) all.times <- log(all.times)
  t_range <- range(all.times)
  y_range <- range(log(do.call(rbind, clldat)$cumhaz))

  par(mfrow = c(1, 1), bty = "n", tcl = -0.15, mgp = c(2.3, 0.5, 0))
  plot(0, 0,
       type = "n", xlab = paste0(ifelse(log_time, "Log-", ""), "Time in ", time_scale),
       ylab = "Log-Cumulative Hazard Rate",
       ylim = y_range, xlim = t_range, yaxt = "n",
       main = paste0(
         "Diagnosis plot for Proportional Hazard assumption\nEndpoint: ", endpoint_name,
         ifelse(subtitle == "", "", "\n"), subtitle
       )
  )
  axis(2, las = 1)

  trts <- names(clldat)
  cols <- c("dodgerblue3", "firebrick3")
  pchs <- c(1, 4)
  for (i in seq_along(clldat)) {
    use_x <- (clldat[[i]]$time / time_unit[[time_scale]])
    if (log_time) use_x <- log(use_x)

    lines(
      y = log(clldat[[i]]$cumhaz),
      x = use_x, col = cols[i]
    )
    points(
      y = log(clldat[[i]]$cumhaz),
      x = use_x,
      col = cols[i], pch = pchs[i], cex = 0.7
    )
  }
  legend("bottomright",
         bty = "n", lty = c(1, 1, 2), cex = 0.8,
         col = cols, pch = pchs, legend = paste0("Treatment: ", trts)
  )
}


#' Plot Schoenfeld residuals for a Cox model fit
#'
#' @param coxobj object returned from \code{\link[survival]{coxph}}
#' @param time_scale a character string, 'year', 'month', 'week' or 'day', time unit of median survival time
#' @param log_time logical, TRUE (default) or FALSE
#' @param endpoint_name a character string, name of the endpoint
#' @param subtitle a character string, subtitle of the plot
#' @examples
#' library(survival)
#' load(system.file("extdata", "combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
#' unweighted_cox <- coxph(Surv(TIME, EVENT == 1) ~ ARM, data = combined_data_tte)
#' resid_plot(unweighted_cox,
#'   time_scale = "month", log_time = TRUE,
#'   endpoint_name = "OS", subtitle = "(Before Matching)"
#' )
#' @return a plot of Schoenfel residuals
#' @export

resid_plot <- function(coxobj, time_scale = "month", log_time = TRUE, endpoint_name = "", subtitle = "") {
  time_unit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)
  if (!time_scale %in% names(time_unit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")

  schresid <- residuals(coxobj, type = "schoenfeld")
  plot_x <- as.numeric(names(schresid)) / time_unit[[time_scale]]
  if (log_time) plot_x <- log(plot_x)
  par(mfrow = c(1, 1), bty = "n", tcl = -0.15, mgp = c(2.3, 0.5, 0))
  plot(schresid ~ plot_x,
       cex = 0.9, col = "navyblue", yaxt = "n",
       ylab = "Unscaled Schoenfeld Residual", xlab = paste0(ifelse(log_time, "Log-", ""), "Time in ", time_scale),
       main = paste0(
         "Diagnostic Plot: Unscaled Schoenfeld Residual\nEndpoint: ", endpoint_name,
         ifelse(subtitle == "", "", "\n"), subtitle
       )
  )
  axis(2, las = 1)
}
