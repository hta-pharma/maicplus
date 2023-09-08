#' Unanchored MAIC for time-to-event endpoint
#'
#' @param useWt a numeric vector of individual MAIC weights, length should the same as \code{nrow(dat)}
#' @param dat a data frame that meet format requirements in 'Details', individual patient data (IPD) of internal trial
#' @param dat_ext a data frame, pseudo IPD from digitized KM curve of external trial
#' @param trt  a character string, name of the interested treatment in internal trial (real IPD)
#' @param trt_ext character string, name of the interested comparator in external trial used to subset \code{dat_ext} (pseudo IPD)
#' @param endpoint_name a character string, name of the endpoint
#' @param time_scale a character string, 'year', 'month', 'week' or 'day', time unit of median survival time
#' @param transform a character string, pass to \code{\link[survival]{cox.zph}}
#'
#' @details Format requirements for input \code{dat} and \code{dat_ext} are to have the following columns
#' \itemize{
#'   \item treatment - character or factor column
#'   \item status - logical column, TRUE for censored/death, FALSE for otherwise
#'   \item time - numeric column, observation time of the \code{status}; unit in days
#' }
#'
#' @return A list of KM plot, analysis table, and diagnostic plot
#' @export

maic_tte_unanchor <- function(useWt, dat, dat_ext, trt, trt_ext,
                              time_scale = "month", endpoint_name = "OS",
                              transform = "log") {
  timeUnit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)

  if (length(useWt) != nrow(dat)) stop("length of useWt should be the same as nrow(dat)")
  if (!time_scale %in% names(timeUnit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")

  res <- list()

  # set up IPD 'dat' with maic weights
  dat$weight <- useWt
  dat <- dat[, c("treatment", "time", "status", "weight")]

  # set up pseudo IPD 'dat_ext' with universal weight of 1
  dat_ext <- dat_ext[dat_ext$treatment == trt_ext, ]
  dat_ext$treatment <- trt_ext
  dat_ext$weight <- 1
  dat_ext <- dat_ext[, names(dat)]

  # merge pseudo IPD and real ipd
  dat <- rbind(dat, dat_ext)
  dat$treatment <- factor(dat$treatment, levels = c(trt_ext, trt))

  # ==> Report 1: KM plot

  # derive km w and w/o weights
  kmobj <- survfit(Surv(time, status) ~ treatment, dat, conf.type = "log-log")
  kmobj_adj <- survfit(Surv(time, status) ~ treatment, dat, weights = dat$weight, conf.type = "log-log")

  par(cex.main = 0.85)
  km_plot(kmobj, kmobj_adj,
    time_scale = time_scale,
    trt = trt, trt_ext = trt_ext,
    endpoint_name = endpoint_name
  )
  res[["plot_km"]] <- grDevices::recordPlot()

  res[["fit_km_data_before"]] <- survfit_makeup(kmobj)
  res[["fit_km_data_after"]] <- survfit_makeup(kmobj_adj)

  # ==> Report 2: Analysis table (Cox model) before and after matching, incl Median Survival Time

  # derive median survival time
  medSurv <- medSurv_makeup(kmobj, legend = "before matching", time_scale = time_scale)
  medSurv_adj <- medSurv_makeup(kmobj_adj, legend = "after matching", time_scale = time_scale)
  medSurv_out <- rbind(medSurv, medSurv_adj)

  res[["report_median_surv"]] <- medSurv_out

  # fit PH Cox regression model
  coxobj <- coxph(Surv(time, status) ~ treatment, dat, robust = T)
  coxobj_adj <- coxph(Surv(time, status) ~ treatment, dat, weights = dat$weight, robust = T)

  res[["fit_cox_model_before"]] <- coxobj
  res[["fit_cox_model_after"]] <- coxobj_adj

  res[["report_overall"]] <- rbind(
    report_table(coxobj, medSurv, tag = paste0("Before/", endpoint_name)),
    report_table(coxobj_adj, medSurv_adj, tag = paste0("After/", endpoint_name))
  )

  # ==> Report 3: Diagnosis Plot

  # grambsch & theaneu ph test
  coxdiag <- cox.zph(coxobj, global = F, transform = transform)
  coxdiag_adj <- cox.zph(coxobj_adj, global = F, transform = transform)

  res[["fit_GT_test_before"]] <- coxdiag
  par(mfrow = c(1, 1), tcl = -0.15)
  plot(coxdiag, yaxt = "n", main = "Grambsch & Terneau Plot (before matching)")
  axis(2, las = 1)
  res[["plot_GT_before"]] <- grDevices::recordPlot()


  res[["fit_GT_test_after"]] <- coxdiag_adj
  par(mfrow = c(1, 1), tcl = -0.15)
  plot(coxdiag_adj, yaxt = "n", main = "Grambsch & Terneau Plot (after matching)")
  axis(2, las = 1)
  res[["plot_GT_after"]] <- grDevices::recordPlot()

  # log-cumulative hazard plot
  log_cum_haz_plot(res[["fit_km_data_before"]],
    time_scale = "month", log_time = TRUE,
    endpoint_name = endpoint_name, subtitle = "(Before Matching)"
  )
  res[["plot_logCH_before"]] <- grDevices::recordPlot()

  log_cum_haz_plot(res[["fit_km_data_after"]],
    time_scale = "month", log_time = TRUE,
    endpoint_name = endpoint_name, subtitle = "(After Matching)"
  )
  res[["plot_logCH_after"]] <- grDevices::recordPlot()

  # schoenfeld residual plot
  resid_plot(coxobj,
    time_scale = "month", log_time = TRUE,
    endpoint_name = endpoint_name, subtitle = "(Before Matching)"
  )
  res[["plot_resid_before"]] <- grDevices::recordPlot()

  resid_plot(coxobj_adj,
    time_scale = "month", log_time = TRUE,
    endpoint_name = endpoint_name, subtitle = "(After Matching)"
  )
  res[["plot_resid_after"]] <- grDevices::recordPlot()

  # output
  res
}

#' helper function: sort out a nice report table to summarize survival analysis results
#' for `maic_tte_unanchor`
#'
#' @param coxobj returned object from \code{\link[survival]{coxph}}
#' @param medSurvobj returned object from \code{\link{medSurv_makeup}}
#' @param tag a string, by default NULL, if specified, an extra 1st column is created in the output
#'
#' @return a data frame with sample size, incidence rate, median survival time with 95% CI, hazard ratio estimate with
#' 95% CI and Wald test of hazard ratio
#'
#' @export

report_table <- function(coxobj, medSurvobj, tag = NULL) {
  hr_res <- format(round(summary(coxobj)$conf.int[-2], 2), nsmall = 2)
  hr_res <- paste0(hr_res[1], "[", hr_res[2], ";", hr_res[3], "]")
  hr_pval <- format(round(summary(coxobj)$waldtest[3], 3), nsmall = 3)
  if (hr_pval == "0.000") hr_pval <- "<0.001"

  meds_report <- format(round(medSurvobj[, c("median", "0.95LCL", "0.95UCL")], 1), nsmall = 1)
  meds_report <- apply(meds_report, 1, function(xx) paste0(xx[1], "[", xx[2], ";", xx[3], "]"))

  desc_res <- cbind(medSurvobj[, "treatment", drop = F],
    data.frame(N = round(medSurvobj$n.max, 1)),
    "n.events(%)" = paste0(
      round(medSurvobj$events, 1), "(",
      format(round(medSurvobj$events * 100 / medSurvobj$n.max, 1), nsmall = 1), ")"
    ),
    "median[95% CI]" = meds_report
  )

  desc_res <- cbind(desc_res[c(2, 1), ], "HR[95% CI]" = c(hr_res, ""), "WaldTest" = c(hr_pval, ""))

  if (!is.null(tag)) desc_res <- cbind(data.frame(Matching = rep(tag, nrow(desc_res))), desc_res)
  desc_res
}
