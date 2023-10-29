#' Anchored MAIC for continuous, binary and time-to-event endpoint
#'
#' This is a wrapper function to provide adjusted effect estimates and relevant statistics in anchored case
#' (i.e. there is a common comparator arm in the internal and external trial).
#'
#' @param ipd_weights an object returned by \code{estimate_weight}
#' @param dat_ipd a data frame that meet format requirements in 'Details', individual patient data (IPD) of internal trial
#' @param ipd_trt_var a string, column name in \code{dat_ipd} that contains the treatment assignment
#' @param dat_pseudo a data frame, pseudo IPD from digitized KM curve of external trial (for time-to-event endpoint) or from contingency table (for binary endpoint)
#' @param pseudo_trt_var a string, column name in \code{dat_ipd} that contains the treatment assignment
#' @param trt_ipd  a string, name of the interested investigation arm in internal trial \code{dat_igd} (real IPD)
#' @param trt_agd a string, name of the interested investigation arm in external trial \code{dat_pseudo} (pseudo IPD)
#' @param trt_common a string, name of the common comparator in internal and external trial
#' @param endpoint_type a string, one out of the following "binary", "tte" (time to event), "continuous"
#' @param eff_measure a string, "Diff" for continuous,"RD" (risk difference), "OR" (odds ratio), "RR" (relative risk) for a binary endpoint; "HR" for a time-to-event endpoint. By default is \code{NULL}, "OR" is used for binary case, otherwise "HR" is used.
#' @param endpoint_name a string, name of time to event endpoint, to be show in the last line of title
#' @param time_scale a string, time unit of median survival time, taking a value of 'year', 'month', 'week' or 'day'
#' @param km_conf_type a string, pass to \code{conf.type} of \code{survfit}
#'
#' @details Format requirements for input \code{dat_ipd} and \code{dat_pseudo} are to have the following columns
#' \itemize{
#'   \item USUBJID - charater, unique subject ID
#'   \item ARM - character or factor column
#'   \item EVENT - numeric, 1 for censored/death, 0 for otherwise
#'   \item TIME - numeric column, observation time of the \code{EVENT}; unit in days
#' }
#'
#' @importFrom grDevices recordPlot
#' @return A list, contains 'descriptive' and 'inferential'
#' @export

maic_anchored <- function(ipd_weights,
                          dat_ipd,
                          ipd_trt_var = "arm",
                          dat_pseudo,
                          pseudo_trt_var = "arm",
                          trt_ipd,
                          trt_agd,
                          trt_common,
                          endpoint_type = "tte",
                          endpoint_name = "Time to Event Endpoint",
                          eff_measure = c("HR","OR","RR","RD","Diff"),
                          # time to event specific args
                          time_scale = "month",
                          km_conf_type = "log-log") {

  # ==> Setup and Precheck ------------------------------------------
  timeUnit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)
  names(dat_ipd) <- toupper(names(dat_ipd))
  names(dat_pseudo) <- toupper(names(dat_pseudo))
  ipd_trt_var <- toupper(ipd_trt_var)
  pseudo_trt_var <- toupper(pseudo_trt_var)
  if (length(eff_measure)>1) eff_measure <- NULL
  if (is.null(eff_measure)) eff_measure <- list(binary="OR",tte="HR",continuous="Diff")[[endpoint_type]]

  # setup ARM column and precheck
  if(!ipd_trt_var %in% names(dat_ipd)) stop("cannot find arm indicator column ipd_trt_var in dat_ipd")
  if(!pseudo_trt_var %in% names(dat_pseudo)) stop("cannot find arm indicator column pseudo_trt_var in dat_pseudo")
  if (ipd_trt_var != "ARM")  dat_ipd$ARM <- dat_ipd[[ipd_trt_var]]
  if (pseudo_trt_var != "ARM") dat_pseudo$ARM <- dat_pseudo[[pseudo_trt_var]]
  dat_ipd$ARM <- as.character(dat_ipd$ARM) # just to avoid potential error when merging
  dat_pseudo$ARM <- as.character(dat_pseudo$ARM) # just to avoid potential error when merging
  if (!trt_ipd %in% dat_ipd$ARM) stop("trt_ipd does not exist in dat_ipd$ARM")
  if (!trt_agd %in% dat_pseudo$ARM) stop("trt_agd does not exist in dat_pseudo$ARM")
  if (!trt_common %in% dat_ipd$ARM) stop("trt_common does not exist in dat_ipd$ARM")
  if (!trt_common %in% dat_pseudo$ARM) stop("trt_common does not exist in dat_pseudo$ARM")
  arms_in_dat_ipd <- unique(dat_ipd$ARM)
  arms_in_dat_pseudo <- unique(dat_ipd$ARM)
  if (length(arms_in_dat_ipd) != 2) stop(paste("In anchored case, there should be two arms in dat_ipd, but you have:", paste(arms_in_dat_ipd, collapse = ",")))
  if (length(arms_in_dat_pseudo) != 2) stop(paste("In anchored case, there should be two arms in dat_pseudo, but you have:", paste(arms_in_dat_pseudo, collapse = ",")))

  # more pre-checks
  endpoint_type <- match.arg(endpoint_type, c("binary","tte","continuous"))
  if(!"maicplus_estimate_weights" %in% class(ipd_weights)) stop("ipd_weights should be an object returned by estimate_weights")
  if (any(duplicated(dat_ipd$USUBJID))) warning("check your dat_ipd, it has duplicated usubjid, this indicates, it might contain multiple endpoints for each subject")
  if (!all(dat_ipd$USUBJID %in% ipd_weights$data$USUBJID)) {
    stop(paste("these pts in dat_ipd cannot be found in ipd_weights",
               paste(setdiff(dat_ipd$USUBJID, ipd_weights$USUBJID), collapse = ","))
    )
  }
  if (!time_scale %in% names(timeUnit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")
  if (endpoint_type == "binary") { # for binary effect measure
    if (any(!c("USUBJID","VALUE") %in% names(dat_ipd))) stop("dat_ipd should have 'USUBJID', 'VALUE' columns at minimum")
    eff_measure <- match.arg(eff_measure, choices = c("OR","RD","RR"), several.ok = FALSE)
  } else if (endpoint_type == "tte"){ # for time to event effect measure
    if(!all(c("USUBJID","TIME", "EVENT", ipd_trt_var) %in% names(dat_ipd))) stop(paste("dat_ipd needs to include at least USUBJID, TIME, EVENT,", ipd_trt_var))
    if(!all(c("TIME", "EVENT", pseudo_trt_var) %in% names(dat_pseudo))) stop(paste("dat_pseudo needs to include at least TIME, EVENT,", pseudo_trt_var))
    eff_measure <- match.arg(eff_measure, choices = c("HR"), several.ok = FALSE)
  } else { # for continuous effect measure
    if (any(!c("USUBJID","VALUE") %in% names(dat_ipd))) stop("dat_ipd should have 'USUBJID', 'VALUE' columns at minimum")
    eff_measure <- match.arg(eff_measure, choices = c("Diff"), several.ok = FALSE)
  }

  # create the hull for the output from this function
  res <- list(
    descriptive = list(),
    inferential = list()
  )

  # prepare ipd and agd data for analysis
  dat_ipd <- dat_ipd[dat_ipd$ARM %in% c(trt_ipd, trt_common), , drop = TRUE]
  dat_pseudo <- dat_pseudo[dat_pseudo$ARM %in% c(trt_agd, trt_common), , drop = TRUE]
  dat_ipd$weights <- ipd_weights$data$weights[match(ipd_weights$data$USUBJID, dat_ipd$USUBJID)]
  dat_pseudo$weights <- 1
  if (!"USUBJID" %in% names(dat_pseudo)) dat_pseudo$USUBJID <- paste0("ID", 1:nrow(dat_pseudo))
  if (any(is.na(dat_ipd$weights))){
    warning(
      paste("these usubjid in dat_ipd have no weight in ipd_weights:",
            paste(dat_ipd$USUBJID[is.na(dat_ipd$weights)], collapse = ","))
    )
  }

  if (endpoint_type == "tte") {
    retain_cols <- c("USUBJID","ARM","TIME","EVENT","weights")
  } else {
    retain_cols <- c("USUBJID","ARM","VALUE","weights")
  }
  dat_ipd <- dat_ipd[, retain_cols, drop = FALSE]
  dat_pseudo <- dat_pseudo[, retain_cols, drop = FALSE]
  dat <- rbind(dat_ipd, dat_pseudo)
  dat_ipd$ARM <- factor(dat_ipd$ARM, levels = c(trt_common, trt_ipd))
  dat_pseudo$ARM <- factor(dat_pseudo$ARM, levels = c(trt_common, trt_agd))
  dat$ARM <- factor(dat$ARM, levels = c(trt_common, trt_agd, trt_ipd))

  # ==> Inferential output ------------------------------------------
  if(endpoint_type == "tte"){
    # Analysis table (Cox model) before and after matching, incl Median Survival Time
    # derive km w and w/o weights
    kmobj_ipd <- survfit(Surv(TIME, EVENT) ~ ARM, dat_ipd, conf.type = km_conf_type)
    kmobj_ipd_adj <- survfit(Surv(TIME, EVENT) ~ ARM, dat_ipd, weights = dat_ipd$weights, conf.type = km_conf_type)
    kmobj_agd <- survfit(Surv(TIME, EVENT) ~ ARM, dat_pseudo, conf.type = km_conf_type)
    res$descriptive[["survfit_ipd_before"]] <- survfit_makeup(kmobj_ipd)
    res$descriptive[["survfit_ipd_after"]] <- survfit_makeup(kmobj_ipd_adj)
    res$descriptive[["survfit_pseudo"]] <- survfit_makeup(kmobj_agd)
    # derive median survival time
    medSurv_ipd <- medSurv_makeup(kmobj_ipd, legend = "IPD, before matching", time_scale = time_scale)
    medSurv_ipd_adj <- medSurv_makeup(kmobj_ipd_adj, legend = "IPD, after matching", time_scale = time_scale)
    medSurv_agd <- medSurv_makeup(kmobj_agd, legend = "AgD, external", time_scale = time_scale)
    medSurv_out <- rbind(medSurv_ipd, medSurv_ipd_adj, medSurv_agd)

    res$inferential[["report_median_surv"]] <- medSurv_out

    # fit PH Cox regression model
    coxobj_ipd <- coxph(Surv(TIME, EVENT) ~ ARM, dat_ipd, robust = T)
    coxobj_ipd_adj <- coxph(Surv(TIME, EVENT) ~ ARM, dat_ipd, weights = weights, robust = T)
    coxobj_agd <- coxph(Surv(TIME, EVENT) ~ ARM, dat_pseudo, robust = T)

    res$inferential[["ipd_coxph_before"]] <- coxobj_ipd
    res$inferential[["ipd_coxph_after"]] <- coxobj_ipd_adj
    res$inferential[["agd_coxph"]] <- coxobj_agd

    # derive ipd exp arm vs agd exp arm via bucher
    res_AC <- as.list(summary(coxobj_ipd_adj)$coef)[c(1,4)]
    res_BC <- as.list(summary(coxobj_agd)$coef)[c(1,4)]
    names(res_AC) <- names(res_BC) <- c("est","se")
    res_AB <- bucher(res_AC, res_BC, conf_lv=0.95)
    res_AB$est <- exp(res_AB$est)
    res_AB$ci_l <- exp(res_AB$ci_l)
    res_AB$ci_u <- exp(res_AB$ci_u)

    res$inferential[["report_overall"]] <- rbind(
      report_table(coxobj_ipd, medSurv_ipd, tag = paste0("IPD/", endpoint_name)),
      report_table(coxobj_ipd_adj, medSurv_ipd_adj, tag = paste0("weighted IPD/", endpoint_name)),
      report_table(coxobj_agd, medSurv_agd, tag = paste0("Agd/", endpoint_name)),
      c(paste0("** adj.",trt_ipd," vs ", trt_agd),
        rep("-",4),
        print_bucher(output = res_AB, pval_digits = 3))
    )
  }

  # output
  res
}
