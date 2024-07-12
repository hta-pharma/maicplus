#' Anchored MAIC for binary and time-to-event endpoint
#'
#' This is a wrapper function to provide adjusted effect estimates and relevant statistics in anchored case (i.e. there
#' is a common comparator arm in the internal and external trial).
#'
#' @param weights_object an object returned by \code{estimate_weight}
#' @param ipd a data frame that meet format requirements in 'Details', individual patient data (IPD) of internal trial
#' @param pseudo_ipd a data frame, pseudo IPD from digitized KM curve of external trial (for time-to-event endpoint) or
#'   from contingency table (for binary endpoint)
#' @param trt_ipd  a string, name of the interested investigation arm in internal trial \code{ipd} (internal IPD)
#' @param trt_agd a string, name of the interested investigation arm in external trial \code{pseudo_ipd} (pseudo IPD)
#' @param trt_common a string, name of the common comparator in internal and external trial
#' @param trt_var_ipd a string, column name in \code{ipd} that contains the treatment assignment
#' @param trt_var_agd a string, column name in \code{ipd} that contains the treatment assignment
#' @param endpoint_type a string, one out of the following "binary", "tte" (time to event)
#' @param endpoint_name a string, name of time to event endpoint, to be show in the last line of title
#' @param time_scale a string, time unit of median survival time, taking a value of 'years', 'months', 'weeks' or
#'   'days'. NOTE: it is assumed that values in TIME column of \code{ipd} and \code{pseudo_ipd} is in the unit of days
#' @param km_conf_type a string, pass to \code{conf.type} of \code{survfit}
#' @param eff_measure a string, "RD" (risk difference), "OR" (odds ratio), "RR" (relative risk) for a binary endpoint;
#'   "HR" for a time-to-event endpoint. By default is \code{NULL}, "OR" is used for binary case, otherwise "HR" is used.
#' @param boot_ci_type a string, one of `c("norm","basic", "stud", "perc", "bca")` to select the type of bootstrap
#'   confidence interval. See [boot::boot.ci] for more details.
#' @param binary_robust_cov_type a string to pass to argument `type` of [sandwich::vcovHC], see possible options in the
#'   documentation of that function. Default is `"HC3"`
#'
#' @details It is required that input \code{ipd} and \code{pseudo_ipd} to have the following
#'   columns. This function is not sensitive to upper or lower case of letters in column names.
#' \itemize{
#'   \item USUBJID - character, unique subject ID
#'   \item ARM - character or factor, treatment indicator, column name does not have to be 'ARM'. User specify in
#'    \code{trt_var_ipd} and \code{trt_var_agd}
#'  }
#'  For time-to-event analysis, the follow columns are required:
#'  \itemize{
#'   \item EVENT - numeric, `1` for censored/death, `0` otherwise
#'   \item TIME - numeric column, observation time of the \code{EVENT}; unit in days
#' }
#' For binary outcomes:
#' \itemize{
#'   \item RESPONSE - numeric, `1` for event occurred, `0` otherwise
#' }
#'
#' @importFrom survival survfit Surv coxph
#' @importFrom lmtest coeftest coefci
#' @importFrom sandwich vcovHC
#' @importFrom boot boot boot.ci
#' @return A list, contains 'descriptive' and 'inferential'
#' @example inst/examples/maic_anchored_binary_ex.R
#' @export

maic_anchored <- function(weights_object,
                          ipd,
                          pseudo_ipd,
                          trt_ipd,
                          trt_agd,
                          trt_common,
                          trt_var_ipd = "ARM",
                          trt_var_agd = "ARM",
                          endpoint_type = "tte",
                          endpoint_name = "Time to Event Endpoint",
                          eff_measure = c("HR", "OR", "RR", "RD"),
                          boot_ci_type = c("norm", "basic", "stud", "perc", "bca"),
                          # time to event specific args
                          time_scale = "months",
                          km_conf_type = "log-log",
                          # binary specific args
                          binary_robust_cov_type = "HC3") {
  # ==> Initial Setup ------------------------------------------
  # ~~~ Create the hull for the output from this function
  res <- list(
    descriptive = list(),
    inferential = list()
  )

  res_AB <- list(
    est = NA,
    se = NA,
    ci_l = NA,
    ci_u = NA,
    pval = NA
  )

  # ~~~ Initial colname process and precheck on effect measure
  names(ipd) <- toupper(names(ipd))
  names(pseudo_ipd) <- toupper(names(pseudo_ipd))
  trt_var_ipd <- toupper(trt_var_ipd)
  trt_var_agd <- toupper(trt_var_agd)
  if (length(eff_measure) > 1) eff_measure <- NULL
  if (is.null(eff_measure)) eff_measure <- list(binary = "OR", tte = "HR")[[endpoint_type]]

  # ~~~ Setup ARM column and make related pre-checks
  if (!trt_var_ipd %in% names(ipd)) stop("cannot find arm indicator column trt_var_ipd in ipd")
  if (!trt_var_agd %in% names(pseudo_ipd)) stop("cannot find arm indicator column trt_var_agd in pseudo_ipd")
  if (trt_var_ipd != "ARM") ipd$ARM <- ipd[[trt_var_ipd]]
  if (trt_var_agd != "ARM") pseudo_ipd$ARM <- pseudo_ipd[[trt_var_agd]]
  ipd$ARM <- as.character(ipd$ARM) # just to avoid potential error when merging

  # ~~~ More pre-checks
  pseudo_ipd$ARM <- as.character(pseudo_ipd$ARM) # just to avoid potential error when merging
  if (!trt_ipd %in% ipd$ARM) stop("trt_ipd does not exist in ipd$ARM")
  if (!trt_agd %in% pseudo_ipd$ARM) stop("trt_agd does not exist in pseudo_ipd$ARM")
  if (!trt_common %in% ipd$ARM) stop("trt_common does not exist in ipd$ARM")
  if (!trt_common %in% pseudo_ipd$ARM) stop("trt_common does not exist in pseudo_ipd$ARM")
  ipd_arms <- unique(ipd$ARM)
  pseudo_ipd_arms <- unique(pseudo_ipd$ARM)
  if (!length(ipd_arms) >= 2) {
    stop("In anchored case, there should be at least two arms in ipd, but you have: ", toString(ipd_arms))
  }
  if (!length(pseudo_ipd_arms) >= 2) {
    stop("In anchored case, there should be at least two arms in pseudo_ipd, but you have: ", toString(pseudo_ipd_arms))
  }
  endpoint_type <- match.arg(endpoint_type, c("binary", "tte"))
  if (!"maicplus_estimate_weights" %in% class(weights_object)) {
    stop("weights_object should be an object returned by estimate_weights")
  }
  if (any(duplicated(ipd$USUBJID))) {
    warning(
      "check your ipd, it has duplicated usubjid, this indicates, ",
      "it might contain multiple endpoints for each subject"
    )
  }
  if (!all(ipd$USUBJID %in% weights_object$data$USUBJID)) {
    stop(
      "These patients in ipd cannot be found in weights_object ",
      toString(setdiff(ipd$USUBJID, weights_object$USUBJID))
    )
  }
  time_scale <- match.arg(arg = time_scale, choices = c("days", "weeks", "months", "years"))
  if (endpoint_type == "binary") { # for binary effect measure

    if (any(!c("USUBJID", "RESPONSE") %in% names(ipd))) stop("ipd should have 'USUBJID', 'RESPONSE' columns at minimum")
    eff_measure <- match.arg(eff_measure, choices = c("OR", "RD", "RR"), several.ok = FALSE)
  } else if (endpoint_type == "tte") { # for time to event effect measure

    if (!all(c("USUBJID", "TIME", "EVENT", trt_var_ipd) %in% names(ipd))) {
      stop("ipd needs to include at least USUBJID, TIME, EVENT, ", toString(trt_var_ipd))
    }
    if (!all(c("TIME", "EVENT", trt_var_agd) %in% names(pseudo_ipd))) {
      stop("pseudo_ipd needs to include at least TIME, EVENT, ", toString(trt_var_agd))
    }
    eff_measure <- match.arg(eff_measure, choices = c("HR"), several.ok = FALSE)
  }
  boot_ci_type <- match.arg(boot_ci_type)

  # ==> IPD and AgD data preparation ------------------------------------------
  # : subset ipd, retain only ipd from interested trts
  ipd <- ipd[ipd$ARM %in% c(trt_ipd, trt_common), , drop = TRUE]
  pseudo_ipd <- pseudo_ipd[pseudo_ipd$ARM %in% c(trt_agd, trt_common), , drop = TRUE]

  # : assign weights to real and pseudo ipd
  ipd$weights <- weights_object$data$weights[match(weights_object$data$USUBJID, ipd$USUBJID)]
  pseudo_ipd$weights <- 1
  if (!"USUBJID" %in% names(pseudo_ipd)) pseudo_ipd$USUBJID <- paste0("ID", seq_len(nrow(pseudo_ipd)))

  # : give warning when individual pts in IPD has no weights
  if (any(is.na(ipd$weights))) {
    ipd <- ipd[!is.na(ipd$weights), , drop = FALSE]
    warning(
      "These USUBJID in IPD have no weight in weights_object and hence excluded from analysis: ",
      toString(ipd$USUBJID[is.na(ipd$weights)])
    )
    if (nrow(ipd) == 0) stop("There are no patients with valid weights in IPD!")
  }

  # : retain necessary columns
  outcome_cols <- if (endpoint_type == "tte") c("TIME", "EVENT") else "RESPONSE"
  retain_cols <- c("USUBJID", "ARM", outcome_cols, "weights")

  ipd <- ipd[, retain_cols, drop = FALSE]
  pseudo_ipd <- pseudo_ipd[, retain_cols, drop = FALSE]

  # : merge real and pseudo ipds, only used if apply contrast method,
  #   since contrast method is not implemented in v0.1, this R obj is not used
  #   just a placeholder
  dat <- rbind(ipd, pseudo_ipd)

  # : setup ARM column as a factor,
  # * these line cannot be move prior to "dat <- rbind(ipd, pseudo_ipd)"
  ipd$ARM <- factor(ipd$ARM, levels = c(trt_common, trt_ipd))
  pseudo_ipd$ARM <- factor(pseudo_ipd$ARM, levels = c(trt_common, trt_agd))
  dat$ARM <- factor(dat$ARM, levels = c(trt_common, trt_agd, trt_ipd))

  # ==> Inferential output ------------------------------------------
  result <- if (endpoint_type == "tte") {
    maic_anchored_tte(
      res,
      res_BC = NULL,
      ipd,
      pseudo_ipd,
      km_conf_type,
      time_scale,
      weights_object,
      endpoint_name,
      trt_ipd,
      trt_agd,
      boot_ci_type
    )
  } else if (endpoint_type == "binary") {
    maic_anchored_binary(
      res,
      res_BC = NULL,
      ipd,
      pseudo_ipd,
      binary_robust_cov_type,
      weights_object,
      endpoint_name,
      eff_measure,
      trt_ipd,
      trt_agd,
      boot_ci_type
    )
  } else {
    stop("Endpoint type ", endpoint_type, " currently unsupported.")
  }

  result
}


# MAIC inference functions for TTE outcome type ------------
maic_anchored_tte <- function(res,
                              res_BC = NULL,
                              ipd,
                              pseudo_ipd,
                              km_conf_type,
                              time_scale,
                              weights_object,
                              endpoint_name,
                              trt_ipd,
                              trt_agd,
                              boot_ci_type) {
  # Analysis table (Cox model) before and after matching, incl Median Survival Time
  # derive km w and w/o weights
  kmobj_ipd <- survfit(Surv(TIME, EVENT) ~ ARM, ipd, conf.type = km_conf_type)
  kmobj_ipd_adj <- survfit(Surv(TIME, EVENT) ~ ARM, ipd, weights = ipd$weights, conf.type = km_conf_type)
  kmobj_agd <- survfit(Surv(TIME, EVENT) ~ ARM, pseudo_ipd, conf.type = km_conf_type)
  res$descriptive[["survfit_ipd_before"]] <- survfit_makeup(kmobj_ipd)
  res$descriptive[["survfit_ipd_after"]] <- survfit_makeup(kmobj_ipd_adj)
  res$descriptive[["survfit_pseudo"]] <- survfit_makeup(kmobj_agd)
  # derive median survival time
  medSurv_ipd <- medSurv_makeup(kmobj_ipd, legend = "IPD, before matching", time_scale = time_scale)
  medSurv_ipd_adj <- medSurv_makeup(kmobj_ipd_adj, legend = "IPD, after matching", time_scale = time_scale)
  medSurv_agd <- medSurv_makeup(kmobj_agd, legend = "AgD, external", time_scale = time_scale)
  medSurv_out <- rbind(medSurv_ipd, medSurv_ipd_adj, medSurv_agd)

  res$inferential[["report_median_surv"]] <- medSurv_out

  # ~~~ Analysis table (Cox model) before and after matching
  # fit PH Cox regression model
  coxobj_ipd <- coxph(Surv(TIME, EVENT) ~ ARM, ipd) # robust = TRUE or not makes a diff
  coxobj_ipd_adj <- coxph(Surv(TIME, EVENT) ~ ARM, ipd, weights = weights, robust = TRUE)
  coxobj_agd <- coxph(Surv(TIME, EVENT) ~ ARM, pseudo_ipd)

  res$inferential[["ipd_model_before"]] <- coxobj_ipd
  res$inferential[["ipd_model_after"]] <- coxobj_ipd_adj
  res$inferential[["agd_model"]] <- coxobj_agd

  # derive ipd exp arm vs agd exp arm via bucher
  res_AC <- as.list(summary(coxobj_ipd_adj)$coef)[c(1, 4)] # est, robust se
  if (is.null(res_BC)) res_BC <- as.list(summary(coxobj_agd)$coef)[c(1, 3)] # est, se
  names(res_AC) <- names(res_BC) <- c("est", "se")

  res_AB <- bucher(res_AC, res_BC, conf_lv = 0.95)
  res_AB$est <- exp(res_AB$est)
  res_AB$ci_l <- exp(res_AB$ci_l)
  res_AB$ci_u <- exp(res_AB$ci_u)

  # : get bootstrapped estimates if applicable
  if (!is.null(weights_object$boot)) {
    keep_rows <- setdiff(seq_len(nrow(weights_object$data)), weights_object$rows_with_missing)
    boot_ipd_id <- weights_object$data[keep_rows, "USUBJID", drop = FALSE]

    boot_ipd <- merge(boot_ipd_id, ipd, by = "USUBJID", all.x = TRUE)
    if (nrow(boot_ipd) != nrow(boot_ipd_id)) stop("ipd has multiple observations for some patients")
    boot_ipd <- boot_ipd[match(boot_ipd$USUBJID, boot_ipd_id$USUBJID), ]

    stat_fun <- function(data, index, w_obj) {
      r <- dynGet("r", ifnotfound = NA) # Get bootstrap iteration
      if (!is.na(r)) {
        if (!all(index == w_obj$boot[, 1, r])) stop("Bootstrap and weight indices don't match")
        boot_ipd <- data[w_obj$boot[, 1, r], ]
        boot_ipd$weights <- w_obj$boot[, 2, r]
      }
      boot_coxobj_dat_adj <- coxph(Surv(TIME, EVENT) ~ ARM, boot_ipd, weights = boot_ipd$weights, robust = TRUE)
      boot_res_AC <- list(est = coef(boot_coxobj_dat_adj)[1], se = sqrt(vcov(boot_coxobj_dat_adj)[1, 1]))
      boot_res_AB <- bucher(boot_res_AC, res_BC, conf_lv = 0.95)
      c(
        est_AB = boot_res_AB$est,
        var_AB = boot_res_AB$se^2,
        se_AB = boot_res_AB$se,
        est_AC = boot_res_AC$est,
        se_AC = boot_res_AC$se,
        var_AC = boot_res_AC$se^2
      )
    }

    # Revert seed to how it was for weight bootstrap sampling
    old_seed <- globalenv()$.Random.seed
    on.exit(suspendInterrupts(set_random_seed(old_seed)))
    set_random_seed(weights_object$boot_seed)

    R <- dim(weights_object$boot)[3]
    boot_res <- boot(boot_ipd, stat_fun, R = R, w_obj = weights_object, strata = weights_object$boot_strata)
    boot_ci <- boot.ci(boot_res, type = boot_ci_type, w_obj = weights_object)

    l_u_index <- switch(boot_ci_type,
      "norm" = list(2, 3, "normal"),
      "basic" = list(4, 5, "basic"),
      "stud" = list(4, 5, "student"),
      "perc" = list(4, 5, "percent"),
      "bca" = list(4, 5, "bca")
    )

    res$inferential[["boot_est"]] <- boot_res
    boot_res_AB <- list(
      est = exp(boot_res$t0[1]),
      se = NA,
      ci_l = exp(boot_ci[[l_u_index[[3]]]][l_u_index[[1]]]),
      ci_u = exp(boot_ci[[l_u_index[[3]]]][l_u_index[[2]]]),
      pval = NA
    )
  } else {
    res$inferential[["boot_est"]] <- NULL
  }

  # : make analysis report table
  res$inferential[["report_overall_robustCI"]] <- rbind(
    report_table_tte(coxobj_ipd, medSurv_ipd, tag = paste0("IPD/", endpoint_name)),
    report_table_tte(coxobj_ipd_adj, medSurv_ipd_adj, tag = paste0("weighted IPD/", endpoint_name)),
    report_table_tte(coxobj_agd, medSurv_agd, tag = paste0("Agd/", endpoint_name)),
    c(
      paste0("** adj.", trt_ipd, " vs ", trt_agd),
      rep("--", 4),
      reformat(res_AB, pval_digits = 3)
    )
  )

  if (is.null(res$inferential[["boot_est"]])) {
    res$inferential[["report_overall_bootCI"]] <- NULL
  } else {
    temp_boot_res <- boot_res_AB
    temp_boot_res$ci_l <- boot_res_AB$ci_l
    temp_boot_res$ci_u <- boot_res_AB$ci_u
    class(temp_boot_res) <- class(res_AB)

    res$inferential[["report_overall_bootCI"]] <- rbind(
      report_table_tte(coxobj_ipd, medSurv_ipd, tag = paste0("IPD/", endpoint_name)),
      report_table_tte(coxobj_ipd_adj, medSurv_ipd_adj, tag = paste0("weighted IPD/", endpoint_name)),
      report_table_tte(coxobj_agd, medSurv_agd, tag = paste0("AgD/", endpoint_name)),
      c(
        paste0("** adj.", trt_ipd, " vs ", trt_agd),
        rep("--", 4),
        reformat(boot_res_AB, pval_digits = 3)
      )
    )
  }

  # output
  res
}

# MAIC inference functions for Binary outcome type ------------
maic_anchored_binary <- function(res,
                                 res_BC = NULL,
                                 ipd,
                                 pseudo_ipd,
                                 binary_robust_cov_type,
                                 weights_object,
                                 endpoint_name,
                                 eff_measure,
                                 trt_ipd,
                                 trt_agd,
                                 boot_ci_type) {
  # ~~~ Analysis table
  # : set up proper link
  glm_link <- switch(eff_measure,
    "RD" = poisson(link = "identity"),
    "RR" = poisson(link = "log"),
    "OR" = binomial(link = "logit")
  )
  res_template <- list(
    est = NA,
    se = NA,
    ci_l = NA,
    ci_u = NA,
    pval = NA
  )

  # : fit glm for binary outcome and robust estimate with weights
  binobj_ipd <- glm(RESPONSE ~ ARM, ipd, family = glm_link)
  binobj_ipd_adj <- glm(RESPONSE ~ ARM, ipd, weights = weights, family = glm_link)
  binobj_agd <- glm(RESPONSE ~ ARM, pseudo_ipd, family = glm_link)

  bin_robust_cov <- sandwich::vcovHC(binobj_ipd_adj, type = binary_robust_cov_type)
  bin_robust_coef <- lmtest::coeftest(binobj_ipd_adj, vcov. = bin_robust_cov)
  bin_robust_ci <- lmtest::coefci(binobj_ipd_adj, vcov. = bin_robust_cov)

  res$inferential[["ipd_model_before"]] <- binobj_ipd
  res$inferential[["ipd_model_after"]] <- binobj_ipd_adj
  res$inferential[["agd_model"]] <- binobj_ipd

  # derive ipd exp arm vs agd exp arm via bucher
  res_AC <- res_template
  res_AC$est <- bin_robust_coef[2, "Estimate"]
  res_AC$se <- bin_robust_coef[2, "Std. Error"]
  res_AC$ci_l <- bin_robust_ci[2, "2.5 %"]
  res_AC$ci_u <- bin_robust_ci[2, "97.5 %"]
  res_AC$pval <- bin_robust_coef[2, "Pr(>|z|)"]

  if (is.null(res_BC)) {
    res_BC <- list(est = coef(binobj_agd)[2], se = sqrt(vcov(binobj_agd)[2, 2]))
  }
  res_AB <- bucher(res_AC, res_BC, conf_lv = 0.95)

  if (eff_measure %in% c("RR", "OR")) {
    res_AB$est <- exp(res_AB$est)
    res_AB$se <- sqrt((exp(res_AB$se^2) - 1) * exp(2 * res_AB$est + res_AB$se^2)) # log normal parameterization
    res_AB$ci_l <- exp(res_AB$ci_l)
    res_AB$ci_u <- exp(res_AB$ci_u)
  } else if (eff_measure == "RD") {
    res_AB$est <- res_AB$est * 100
    res_AB$se <- res_AB$se * 100
    res_AB$ci_l <- res_AB$ci_l * 100
    res_AB$ci_u <- res_AB$ci_u * 100
  }

  # : get bootstrapped estimates if applicable
  if (!is.null(weights_object$boot)) {
    keep_rows <- setdiff(seq_len(nrow(weights_object$data)), weights_object$rows_with_missing)
    boot_ipd_id <- weights_object$data[keep_rows, "USUBJID", drop = FALSE]

    boot_ipd <- merge(boot_ipd_id, ipd, by = "USUBJID", all.x = TRUE)
    if (nrow(boot_ipd) != nrow(boot_ipd_id)) stop("ipd has multiple observations for some patients")
    boot_ipd <- boot_ipd[match(boot_ipd$USUBJID, boot_ipd_id$USUBJID), ]

    stat_fun <- function(data, index, w_obj, eff_measure) {
      r <- dynGet("r", ifnotfound = NA) # Get bootstrap iteration
      if (!is.na(r)) {
        if (!all(index == w_obj$boot[, 1, r])) stop("Bootstrap and weight indices don't match")
        boot_ipd <- data[w_obj$boot[, 1, r], ]
        boot_ipd$weights <- w_obj$boot[, 2, r]
      }

      boot_binobj_dat_adj <- glm(RESPONSE ~ ARM, boot_ipd, weights = boot_ipd$weights, family = glm_link)
      boot_AC_est <- coef(boot_binobj_dat_adj)[2]
      boot_AC_var <- vcov(boot_binobj_dat_adj)[2, 2]

      boot_res_AC <- list(est = boot_AC_est, se = sqrt(boot_AC_var))
      boot_res_AB <- bucher(boot_res_AC, res_BC, conf_lv = 0.95)

      c(
        est_AB = boot_res_AB$est,
        var_AB = boot_res_AB$se^2,
        se_AB = boot_res_AB$se,
        est_AC = boot_res_AC$est,
        se_AC = boot_res_AC$se,
        var_AC = boot_res_AC$se^2
      )
    }

    # Revert seed to how it was for weight bootstrap sampling
    old_seed <- globalenv()$.Random.seed
    on.exit(suspendInterrupts(set_random_seed(old_seed)))
    set_random_seed(weights_object$boot_seed)

    R <- dim(weights_object$boot)[3]
    boot_res <- boot(
      boot_ipd,
      stat_fun,
      R = R,
      w_obj = weights_object,
      eff_measure = eff_measure,
      strata = weights_object$boot_strata
    )
    boot_ci <- boot.ci(boot_res, type = boot_ci_type, w_obj = weights_object)

    l_u_index <- switch(boot_ci_type,
      "norm" = list(2, 3, "normal"),
      "basic" = list(4, 5, "basic"),
      "stud" = list(4, 5, "student"),
      "perc" = list(4, 5, "percent"),
      "bca" = list(4, 5, "bca")
    )

    transform_estimate <- switch(eff_measure,
      "RD" = function(x) x * 100,
      "RR" = exp,
      "OR" = exp
    )

    res$inferential[["boot_est"]] <- boot_res
    boot_res_AB <- list(
      est = transform_estimate(boot_res$t0[1]),
      se = NA,
      ci_l = transform_estimate(boot_ci[[l_u_index[[3]]]][l_u_index[[1]]]),
      ci_u = transform_estimate(boot_ci[[l_u_index[[3]]]][l_u_index[[2]]]),
      pval = NA
    )
  } else {
    res$inferential[["boot_est"]] <- NULL
  }

  # : make analysis report table
  tags <- paste0(c("IPD/", "weighted IPD/", "AgD/"), endpoint_name)
  res$inferential[["report_overall_robustCI"]] <- rbind(
    report_table_binary(binobj_ipd, tag = tags[1], eff_measure = eff_measure),
    report_table_binary(binobj_ipd_adj, res_AB, tag = tags[2], eff_measure = eff_measure),
    report_table_binary(binobj_agd, tag = tags[3], eff_measure = eff_measure),
    c(
      paste0("** adj.", trt_ipd, " vs ", trt_agd),
      rep("--", 3),
      reformat(res_AB, pval_digits = 3)
    )
  )

  if (is.null(res$inferential[["boot_est"]])) {
    res$inferential[["report_overall_bootCI"]] <- NULL
  } else {
    res$inferential[["report_overall_bootCI"]] <- rbind(
      report_table_binary(binobj_ipd, tag = tags[1], eff_measure = eff_measure),
      report_table_binary(binobj_ipd_adj, boot_res_AB, tag = tags[2], eff_measure = eff_measure),
      report_table_binary(binobj_agd, tag = tags[3], eff_measure = eff_measure),
      c(
        paste0("** adj.", trt_ipd, " vs ", trt_agd),
        rep("--", 3),
        reformat(boot_res_AB, pval_digits = 3)
      )
    )
  }

  # output
  res
}
