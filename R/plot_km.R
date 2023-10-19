#' Kaplan-Meier plot of before and after adjustments
#'
#' @param km_fit_before returned object from \code{survival::survfit} before adjustment
#' @param km_fit_after returned object from \code{survival::survfit} after adjustment
#' @param time_scale time unit of median survival time, taking a value of 'year', 'month', 'week' or 'day'
#' @param trt internal trial treatment
#' @param trt_ext external trial treatment
#' @param endpoint_name name of the endpoint
#' @param line_col color of the line curves with the order of external, internal unadjusted, and internal adjusted
#'
#' @return a Kaplan-Meier plot of before and after adjustments
#' @examples
#' library(survival)
#' load(system.file("extdata", "combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
#' kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte, conf.type = "log-log")
#' kmobj_adj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte,
#'   weights = combined_data_tte$weights, conf.type = "log-log"
#' )
#' par(cex.main = 0.85)
#' km_plot(kmobj, kmobj_adj, time_scale = "month", trt = "A", trt_ext = "B", endpoint_name = "OS")
#' @export

km_plot <- function(km_fit_before, km_fit_after = NULL, time_scale, trt, trt_ext, endpoint_name = "",
                    line_col = c("#5450E4", "#00857C", "#6ECEB2")) {
  time_unit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)

  if (!time_scale %in% names(time_unit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")

  # prepare data for plot
  pd_be <- survfit_makeup(km_fit_before)
  if (!is.null(km_fit_after)) pd_af <- survfit_makeup(km_fit_after)

  # set up x axis (time)
  if (!is.null(km_fit_after)) {
    max_t <- max(km_fit_before$time, km_fit_after$time)
  } else {
    max_t <- max(km_fit_before$time)
  }
  t_range <- c(0, (max_t / time_unit[[time_scale]]) * 1.07)

  # base plot
  par(mfrow = c(1, 1), bty = "n", tcl = -0.15, mgp = c(2.3, 0.5, 0))
  plot(0, 0,
       type = "n", xlab = paste0("Time in ", time_scale), ylab = "Survival Probability",
       ylim = c(0, 1), xlim = t_range, yaxt = "n",
       main = paste0(
         "Kaplan-Meier Curves of Comparator ", ifelse(!is.null(km_fit_after), "(AgD) ", ""),
         "and Treatment ", ifelse(!is.null(km_fit_after), "(IPD)", ""),
         "\nEndpoint: ", endpoint_name
       )
  )
  axis(2, las = 1)

  # add km lines from external trial
  lines(
    y = pd_be[[trt_ext]]$surv,
    x = (pd_be[[trt_ext]]$time / time_unit[[time_scale]]), col = line_col[1],
    type = "s"
  )
  tmpid <- pd_be[[trt_ext]]$censor == 1
  points(
    y = pd_be[[trt_ext]]$surv[tmpid],
    x = (pd_be[[trt_ext]]$time[tmpid] / time_unit[[time_scale]]),
    col = line_col[1], pch = 3, cex = 0.7
  )

  # add km lines from internal trial before adjustment
  lines(
    y = pd_be[[trt]]$surv,
    x = (pd_be[[trt]]$time / time_unit[[time_scale]]), col = line_col[2],
    type = "s"
  )
  tmpid <- pd_be[[trt]]$censor == 1
  points(
    y = pd_be[[trt]]$surv[tmpid],
    x = (pd_be[[trt]]$time[tmpid] / time_unit[[time_scale]]),
    col = line_col[2], pch = 3, cex = 0.7
  )

  # add km lines from internal trial after adjustment
  if (!is.null(km_fit_after)) {
    lines(
      y = pd_af[[trt]]$surv,
      x = (pd_af[[trt]]$time / time_unit[[time_scale]]), col = line_col[3], lty = 2,
      type = "s"
    )
    tmpid <- pd_af[[trt]]$censor == 1
    points(
      y = pd_af[[trt]]$surv[tmpid],
      x = (pd_af[[trt]]$time[tmpid] / time_unit[[time_scale]]),
      col = line_col[3], pch = 3, cex = 0.7
    )
  }

  use_leg <- 1:ifelse(is.null(km_fit_after), 2, 3)
  # add legend
  legend("topright",
         bty = "n", lty = c(1, 1, 2)[use_leg], cex = 0.8, col = line_col[use_leg],
         legend = c(
           paste0("Comparator: ", trt_ext),
           paste0("Treatment: ", trt),
           paste0("Treatment: ", trt, " (with weights)")
         )[use_leg]
  )
}
