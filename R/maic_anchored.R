#' Anchored MAIC for continuous, binary and time-to-event endpoint
#'
#' This is a wrapper function to provide adjusted effect estimates and relevant statistics in anchored case
#' (i.e. there is a common comparator arm in the internal and external trial).
#'
#' @param ipd_weights a numeric vector of individual MAIC weights, length should the same as \code{nrow(dat_igd)}
#' @param tte_dat_ipd a data frame that meet format requirements in 'Details', individual patient data (IPD) of internal trial
#' @param ipd_trt_var a string, column name in \code{dat_ipd} that contains the treatment assignment
#' @param tte_dat_pseudo a data frame, pseudo IPD from digitized KM curve of external trial (for time-to-event endpoint) or from contingency table (for binary endpoint)
#' @param pseudo_trt_var a string, column name in \code{dat_ipd} that contains the treatment assignment
#' @param trt_ipd  a string, name of the interested investigational arm in internal trial \code{dat_igd} (real IPD)
#' @param trt_agd a string, name of the interested investigational arm in external trial \code{dat_pseudo} (pseudo IPD)
#' @param trt_common a string, name of the common comparator in internal and external trial
#' @param endpoint_type a string, one out of the following "binary", "tte" (time to event), "continuous"
#' @param endpoint_name a string, name of the endpoint, used in the reported tables and graphs
#' @param eff_measure a string, "Diff" for continuous,"RD" (risk difference), "OR" (odds ratio), "RR" (relative risk) for a binary endpoint; "HR" for a time-to-event endpoint. By default is \code{NULL}, "OR" is used for binary case, otherwise "HR" is used.
#' @param time_scale a string, 'year', 'month', 'week' or 'day', time unit of median survival time; only relevant when \code{binary}=FALSE
#' @param transform a string, pass to \code{\link[survival]{cox.zph}}; only relevant when \code{endpoint_type}='tte'
#' @param full_output logical, if FALSE, KM plot and PH diagnosis (for time-to-event endpoint) will be muted. The function will return only a simplified table of adjusted effect estimates; only relevant when \code{endpoint_type}='tte'
#'
#' @details Format requirements for input \code{dat_ipd} and \code{dat_pseudo} are to have the following columns
#' \itemize{
#'   \item arm - character or factor column
#'   \item status - logical column, TRUE for censored/death, FALSE for otherwise
#'   \item time - numeric column, observation time of the \code{status}; unit in days
#' }
#'
#' @importFrom grDevices recordPlot
#' @return A list of KM plot, analysis table, and diagnostic plot
#' @export

maic_anchored <- function(ipd_weights,
                          tte_dat_ipd,
                          ipd_trt_var = "arm",
                          tte_dat_pseudo,
                          pseudo_trt_var = "arm",
                          trt_ipd,
                          trt_agd,
                          trt_common,
                          endpoint_type = "tte",
                          endpoint_name = "OS",
                          eff_measure = c("HR","OR","RR","RD","Diff"),
                          full_output = TRUE,
                          # time to event specific args
                          time_scale = "month",
                          km_plot_type = c("basic","ggplot"),
                          km_conf_type = "log-log",
                          km_layout = c("all","by_trial","by_arm"),
                          transform = "log"
                          ) {

  timeUnit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)

  # pre-check
  endpoint_type <- match.arg(endpoint_type, c("binary","tte","continuous"))
  km_layout <- match.arg(km_layout, c("all","by_trial","by_arm"))
  km_plot_type <- match.arg(km_plot_type, c("basic","ggplot"))
  if (length(useWt) != nrow(dat_ipd)) stop("length of useWt should be the same as nrow(dat)")
  if (!time_scale %in% names(timeUnit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")
  if (endpoint_type == "binary") {
    if (any(!c("value", "weight") %in% names(dat_ipd))) stop("dat_ipd should have 'value','weight' columns at minimum")
    eff_measure <- match.arg(eff_measure, choices = c("OR","RD","RR"), several.ok = FALSE)
  } else if (endpoint_type == "tte"){
    if (any(!c("time", "status", "weight") %in% names(dat_ipd))) stop("dat_ipd should have 'time','status','weight' columns at minimum")
    eff_measure <- match.arg(eff_measure, choices = c("HR"), several.ok = FALSE)
  } else {
    if (any(!c("value", "weight") %in% names(dat_ipd))) stop("dat_ipd should have 'value','weight' columns at minimum")
    eff_measure <- match.arg(eff_measure, choices = c("Diff"), several.ok = FALSE)
  }
  if ("usubjid" %in% names(dat_ipd)) stop("dat_ipd should contain USUBJID column, it is used to find the right weights from useWt")
  if ("usubjid" %in% names(useWt)) stop("useWt should contain USUBJID column, it is used to find the right weights for dat_ipd")
  if (any(duplicated(dat_ipd$usubjid))) stop("check your dat_ipd, it has duplicated usubjid, this indicates, it might contain multiple endpoints for each subject")
  if (!trt_ipd %in% dat_ipd$arm) stop("trt_ipd does not exist in dat_ipd$arm")
  if (!trt_agd %in% dat_pseudo$arm) stop("trt_agd does not exist in dat_pseudo$arm")
  if (!trt_common %in% dat_ipd$arm) stop("trt_common does not exist in dat_ipd$arm")
  if (!trt_common %in% dat_pseudo$arm) stop("trt_common does not exist in dat_pseudo$arm")
  arms_in_dat_ipd <- unique(dat_ipd$arm)
  arms_in_dat_pseudo <- unique(dat_ipd$arm)
  if (length(arms_in_dat_ipd) != 2) stop(paste("In anchored case, there should be two arms in dat_ipd, but you have:", paste(arms_in_dat_ipd, collapse = ",")))
  if (length(arms_in_dat_pseudo) != 2) stop(paste("In anchored case, there should be two arms in dat_pseudo, but you have:", paste(arms_in_dat_pseudo, collapse = ",")))

  # setup
  if (length(eff_measure)>1) eff_measure <- NULL
  if (is.null(eff_measure)) eff_measure <- list(binary="OR",tte="HR",continuous="Diff")[[endpoint_type]]
  km_layout <- match.arg(km_layout)
  if(!ipd_trt_var %in% names(dat_ipd)) stop("cannot find arm indicator column ipd_trt_var in dat_ipd")
  if(!pseudo_trt_var %in% names(dat_agd)) stop("cannot find arm indicator column pseudo_trt_var in dat_pseudo")
  dat_ipd$arm <- dat_ipd[[ipd_trt_var]]
  dat_pseudo$arm <- dat_pseudo[[pseudo_trt_var]]
  names(dat_ipd) <- tolower(names(dat_ipd))
  names(dat_pseudo) <- tolower(names(dat_pseudo))
  names(useWt) <- tolower(names(useWt))

  # create the hull for the output from this function
  res <- list(
    descriptive = list(),
    inferential = list(),
    diagnostics = list()
  )

  # set up IPD 'dat' with maic weights
  dat_ipd$weight <- useWt$weights[match(useWt$usubjid, dat_ipd$usubjid),,drop=FALSE]
  if (any(is.na(dat_ipd$weight))){
    warning(
      paste("these usubjid in dat_ipd have no weight in useWt:",
            paste(dat_ipd$usubjid[is.na(dat_ipd$weight)], collapse = ","))
    )
  }
  if (binary) {
    dat_ipd <- dat_ipd[, c("arm", "value", "weight")]
  } else {
    dat_ipd <- dat_ipd[, c("arm", "time", "status", "weight")]
  }
  dat_ipd <- dat_ipd[dat_ipd$arm %in% c(trt_ipd, trt_common), , drop = FALSE]
  dat_ipd$arm <- as.character(dat_ipd$arm) # just to avoid potential error when merging

  # set up pseudo IPD 'dat_pseudo' with universal weight of 1
  dat_pseudo <- dat_pseudo[dat_pseudo$arm %in% c(trt_agd, trt_common), , drop = FALSE]
  dat_pseudo$weight <- 1
  if (binary) {
    dat_pseudo <- dat_pseudo[, c("arm", "value", "weight")]
  } else {
    dat_pseudo <- dat_pseudo[, c("arm", "time", "status", "weight")]
  }
  dat_pseudo$arm <- as.character(dat_pseudo$arm) # just to avoid potential error when merging

  # merge pseudo IPD and real ipd
  dat <- rbind(dat_ipd, dat_pseudo)
  dat$arm <- factor(dat$arm, levels = c(trt_common, trt_agd, trt_ipd))

  # ==> Descriptive output ------------------------------------------
  if(full_output){
    # ** KM plot
    # derive km w and w/o weights
    kmobj_ipd <- survfit(Surv(time, status) ~ arm, dat_ipd, conf.type = km_conf_type)
    kmobj_ipd_adj <- survfit(Surv(time, status) ~ arm, dat_ipd, weights = dat_ipd$weight, conf.type = km_conf_type)
    kmobj_pseudo <- survfit(Surv(time, status) ~ arm, dat_pseudo, conf.type = km_conf_type)
    # plotting
    par(cex.main = 0.85)
    km_plot(kmobj, kmobj_adj,
            time_scale = time_scale,
            trt = trt, trt_ext = trt_ext,
            endpoint_name = endpoint_name
    )
    # save result
    res$descriptive[["plot_km"]] <- grDevices::recordPlot()
    res$descriptive[["survfit_ipd_before"]] <- survfit_makeup(kmobj_ipd)
    res$descriptive[["survfit_ipd_after"]] <- survfit_makeup(kmobj_ipd_adj)
    res$descriptive[["survfit_pseudo"]] <- survfit_makeup(kmobj_pseudo)
  }


  # ==> Inferential output ------------------------------------------
  # ==> Report 2: Analysis table (Cox model) before and after matching, incl Median Survival Time

  # derive median survival time
  medSurv <- medSurv_makeup(kmobj, legend = "before matching", time_scale = time_scale)
  medSurv_adj <- medSurv_makeup(kmobj_adj, legend = "after matching", time_scale = time_scale)
  medSurv_out <- rbind(medSurv, medSurv_adj)

  res[["report_median_surv"]] <- medSurv_out

  # fit PH Cox regression model
  coxobj <- coxph(Surv(time, status) ~ arm, dat, robust = T)
  coxobj_adj <- coxph(Surv(time, status) ~ arm, dat, weights = dat$weight, robust = T)

  res[["fit_cox_model_before"]] <- coxobj
  res[["fit_cox_model_after"]] <- coxobj_adj

  res[["report_overall"]] <- rbind(
    report_table(coxobj, medSurv, tag = paste0("Before/", endpoint_name)),
    report_table(coxobj_adj, medSurv_adj, tag = paste0("After/", endpoint_name))
  )

  # ==> Diagnosis output ------------------------------------------
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
