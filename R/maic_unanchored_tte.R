#' Unanchored MAIC for time-to-event endpoint
#'
#' @param useWt a numeric vector of individual MAIC weights, length should the same as \code{nrow(dat)}
#' @param dat a data frame that meet format requirements in 'Details', individual patient data (IPD) of internal trial
#' @param dat_ext a data frame, pseudo IPD from digitized KM curve of external trial
#' @param trt  a character string, name of the interested treatment in internal trial (real IPD)
#' @param trt_ext character string, name of the interested comparator in external trial used to subset
#' \code{dat_ext} (pseudo IPD)
#' @param endpoint_name a character string, name of the endpoint
#' @param time_scale a character string, 'years', 'months', 'weeks' or 'days', time unit of median survival time
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
#' @importFrom grDevices recordPlot
#' @export

maic_tte_unanchor <- function(useWt, dat, dat_ext, trt, trt_ext,
                              time_scale = "months", endpoint_name = "OS",
                              transform = "log") {
  if (length(useWt) != nrow(dat)) stop("length of useWt should be the same as nrow(dat)")
  time_scale <- match.arg(arg = time_scale, choices = c("days", "weeks", "months", "years"))

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
  coxobj <- coxph(Surv(time, status) ~ treatment, dat, robust = TRUE)
  coxobj_adj <- coxph(Surv(time, status) ~ treatment, dat, weights = dat$weight, robust = TRUE)

  res[["fit_cox_model_before"]] <- coxobj
  res[["fit_cox_model_after"]] <- coxobj_adj

  res[["report_overall"]] <- rbind(
    report_table(coxobj, medSurv, tag = paste0("Before/", endpoint_name)),
    report_table(coxobj_adj, medSurv_adj, tag = paste0("After/", endpoint_name))
  )

  # ==> Report 3: Diagnosis Plot

  # grambsch & theaneu ph test
  coxdiag <- cox.zph(coxobj, global = FALSE, transform = transform)
  coxdiag_adj <- cox.zph(coxobj_adj, global = FALSE, transform = transform)

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
  ph_diagplot_lch(res[["fit_km_data_before"]],
    time_scale = time_scale, log_time = TRUE,
    endpoint_name = endpoint_name, subtitle = "(Before Matching)"
  )
  res[["plot_logCH_before"]] <- grDevices::recordPlot()

  ph_diagplot_lch(res[["fit_km_data_after"]],
    time_scale = time_scale, log_time = TRUE,
    endpoint_name = endpoint_name, subtitle = "(After Matching)"
  )
  res[["plot_logCH_after"]] <- grDevices::recordPlot()

  # schoenfeld residual plot
  ph_diagplot_schoenfeld(coxobj,
    time_scale = time_scale, log_time = TRUE,
    endpoint_name = endpoint_name, subtitle = "(Before Matching)"
  )
  res[["plot_resid_before"]] <- grDevices::recordPlot()

  ph_diagplot_schoenfeld(coxobj_adj,
    time_scale = time_scale, log_time = TRUE,
    endpoint_name = endpoint_name, subtitle = "(After Matching)"
  )
  res[["plot_resid_after"]] <- grDevices::recordPlot()

  # output
  res
}

