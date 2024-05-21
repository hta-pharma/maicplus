#' Unanchored MAIC for binary and time-to-event endpoint
#'
#' This is a wrapper function to provide adjusted effect estimates and relevant statistics in unanchored case (i.e.
#' there is no common comparator arm in the internal and external trial).
#'
#' @param weights_object an object returned by \code{estimate_weight}
#' @param ipd a data frame that meet format requirements in 'Details', individual patient data (IPD) of internal trial
#' @param pseudo_ipd a data frame, pseudo IPD from digitized KM curve of external trial (for time-to-event endpoint) or
#'   from contingency table (for binary endpoint)
#' @param trt_ipd  a string, name of the interested investigation arm in internal trial \code{dat_igd} (real IPD)
#' @param trt_agd a string, name of the interested investigation arm in external trial \code{pseudo_ipd} (pseudo IPD)
#' @param trt_var_ipd a string, column name in \code{ipd} that contains the treatment assignment
#' @param trt_var_agd a string, column name in \code{ipd} that contains the treatment assignment
#' @param endpoint_type a string, one out of the following "binary", "tte" (time to event)
#' @param eff_measure a string, "RD" (risk difference), "OR" (odds ratio), "RR" (relative risk) for a binary endpoint;
#'   "HR" for a time-to-event endpoint. By default is \code{NULL}, "OR" is used for binary case, otherwise "HR" is used.
#' @param boot_ci_is_quantile a logical, specify if the 95% bootstrapped confidence interval should be derived by sample
#'   quantile. Default FALSE, which the estimates assumes to follow asymptotic normal (only if `eff_measure` is "RD") or
#'   log-normal with a variance that can be approximated by bootstrapped sample of the estimate. This default option may
#'   be handy when the number of bootstrap iterations is not big.
#' @param endpoint_name a string, name of time to event endpoint, to be show in the last line of title
#' @param time_scale a string, time unit of median survival time, taking a value of 'years', 'months', 'weeks' or
#'   'days'. NOTE: it is assumed that values in TIME column of \code{ipd} and \code{pseudo_ipd} is in the unit of days
#' @param km_conf_type a string, pass to \code{conf.type} of \code{survfit}
#' @param binary_robust_cov_type a string to pass to argument "type" of \link[clubSandwich]{vcovCR}, see viable options
#'   in the documentation of that function. Default is "CR2"
#'
#' @details For time-to-event analysis, it is required that input \code{ipd} and \code{pseudo_ipd} to have the following
#'   columns. This function is not sensitive to upper or lower case of letters in column names.
#' \itemize{
#'   \item USUBJID - character, unique subject ID
#'   \item ARM - character or factor, treatment indicator, column name does not have to be 'ARM'. User specify in
#'   \code{trt_var_ipd} and \code{trt_var_agd}
#'   \item EVENT - numeric, 1 for censored/death, 0 for otherwise
#'   \item TIME - numeric column, observation time of the \code{EVENT}; unit in days
#' }
#'
#' @importFrom survival survfit Surv
#' @return A list, contains 'descriptive' and 'inferential'
#' @export

maic_unanchored <- function(weights_object,
                            ipd,
                            pseudo_ipd,
                            trt_ipd,
                            trt_agd,
                            trt_var_ipd = "ARM",
                            trt_var_agd = "ARM",
                            endpoint_type = "tte",
                            endpoint_name = "Time to Event Endpoint",
                            eff_measure = c("HR", "OR", "RR", "RD"),
                            boot_ci_is_quantile = FALSE,
                            # time to event specific args
                            time_scale = "months",
                            km_conf_type = "log-log",
                            # binary specific args
                            binary_robust_cov_type = "CR2") {
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
  pseudo_ipd$ARM <- as.character(pseudo_ipd$ARM) # just to avoid potential error when merging
  if (!trt_ipd %in% ipd$ARM) stop("trt_ipd does not exist in ipd$ARM")
  if (!trt_agd %in% pseudo_ipd$ARM) stop("trt_agd does not exist in pseudo_ipd$ARM")

  # ~~~ More pre-checks
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
      "These pts in ipd cannot be found in weights_object ",
      toString(setdiff(ipd$USUBJID, weights_object$USUBJID))
    )
  }
  time_scale <- match.arg(arg = time_scale, choices = c("days", "weeks", "months", "years"))
  if (endpoint_type == "binary") { # for binary effect measure

    if (any(!c("USUBJID", "RESPONSE") %in% names(ipd))) stop("ipd should have 'USUBJID', 'RESPONSE' columns at minimum")
    eff_measure <- match.arg(eff_measure, choices = c("OR", "RD", "RR"), several.ok = FALSE)
    # !! add a check of binary_robust_cov_type
  } else if (endpoint_type == "tte") { # for time to event effect measure

    if (!all(c("USUBJID", "TIME", "EVENT", trt_var_ipd) %in% names(ipd))) {
      stop("ipd needs to include at least USUBJID, TIME, EVENT, ", trt_var_ipd)
    }
    if (!all(c("TIME", "EVENT", trt_var_agd) %in% names(pseudo_ipd))) {
      stop("pseudo_ipd needs to include at least TIME, EVENT, ", trt_var_agd)
    }
    eff_measure <- match.arg(eff_measure, choices = c("HR"), several.ok = FALSE)
  }

  # ==> IPD and AgD data preparation ------------------------------------------
  # : subset ipd, retain only ipd from interested trts
  ipd <- ipd[ipd$ARM == trt_ipd, , drop = TRUE]
  pseudo_ipd <- pseudo_ipd[pseudo_ipd$ARM == trt_agd, , drop = TRUE]

  # : assign weights to real and pseudo ipd
  ipd$weights <- weights_object$data$weights[match(weights_object$data$USUBJID, ipd$USUBJID)]
  pseudo_ipd$weights <- 1

  # : necessary formatting for pseudo ipd
  if (!"USUBJID" %in% names(pseudo_ipd)) pseudo_ipd$USUBJID <- paste0("ID", seq_len(nrow(pseudo_ipd)))
  if ("RESPONSE" %in% names(pseudo_ipd) && is.logical(pseudo_ipd$RESPONSE)) {
    pseudo_ipd$RESPONSE <- as.numeric(pseudo_ipd$RESPONSE)
  }

  # : give warning when individual pts in IPD has no weights
  if (any(is.na(ipd$weights))) {
    ipd <- ipd[!is.na(ipd$weights), , drop = FALSE]
    warning(
      paste(
        "these usubjid in ipd have no weight in weights_object, and hence excluded from analysis:",
        paste(ipd$USUBJID[is.na(ipd$weights)], collapse = ",")
      )
    )
    if (nrow(ipd) == 0) stop("there is no pts with weight in IPD!!")
  }

  # : retain necessary columns
  if (endpoint_type == "tte") {
    retain_cols <- c("USUBJID", "ARM", "TIME", "EVENT", "weights")
  } else {
    retain_cols <- c("USUBJID", "ARM", "RESPONSE", "weights")
  }
  ipd <- ipd[, retain_cols, drop = FALSE]
  pseudo_ipd <- pseudo_ipd[, retain_cols, drop = FALSE]

  # : merge real and pseudo ipds
  dat <- rbind(ipd, pseudo_ipd)
  dat$ARM <- factor(dat$ARM, levels = c(trt_agd, trt_ipd))

  # ==> Inferential output ------------------------------------------

  result <- if (endpoint_type == "tte") {
    maic_unanchored_tte(
      res, res_AB, dat, ipd, pseudo_ipd, km_conf_type, time_scale,
      weights_object, endpoint_name, boot_ci_is_quantile, trt_ipd, trt_agd
    )
  } else if (endpoint_type == "binary") {
    maic_unanchored_binary(
      res, res_AB, dat, ipd, pseudo_ipd, binary_robust_cov_type,
      weights_object, endpoint_name, eff_measure, boot_ci_is_quantile, trt_ipd, trt_agd
    )
  } else {
    stop("Endpoint type ", endpoint_type, " currently unsupported.")
  }

  # output
  result
}

# MAIC inference functions for TTE outcome type ------------

maic_unanchored_tte <- function(res,
                                res_AB,
                                dat,
                                ipd,
                                pseudo_ipd,
                                km_conf_type,
                                time_scale,
                                weights_object,
                                endpoint_name,
                                boot_ci_is_quantile,
                                trt_ipd,
                                trt_agd) {
  # ~~~ Descriptive table before and after matching
  # : derive km w and w/o weights
  kmobj_dat <- survfit(Surv(TIME, EVENT) ~ ARM, dat, conf.type = km_conf_type)
  kmobj_dat_adj <- survfit(Surv(TIME, EVENT) ~ ARM, dat, weights = dat$weights, conf.type = km_conf_type)
  res$descriptive[["survfit_before"]] <- survfit_makeup(kmobj_dat)
  res$descriptive[["survfit_after"]] <- survfit_makeup(kmobj_dat_adj)
  # : derive median survival time
  medSurv_dat <- medSurv_makeup(kmobj_dat, legend = "Before matching", time_scale = time_scale)
  medSurv_dat_adj <- medSurv_makeup(kmobj_dat_adj, legend = "After matching", time_scale = time_scale)
  medSurv_out <- rbind(medSurv_dat, medSurv_dat_adj)

  res$inferential[["report_median_surv"]] <- medSurv_out

  # ~~~ Analysis table (Cox model) before and after matching
  # : fit PH Cox regression model
  coxobj_dat <- coxph(Surv(TIME, EVENT) ~ ARM, dat)
  coxobj_dat_adj <- coxph(Surv(TIME, EVENT) ~ ARM, dat, weights = weights, robust = TRUE)

  res$inferential[["coxph_before"]] <- coxobj_dat
  res$inferential[["coxph_after"]] <- coxobj_dat_adj

  # : derive ipd exp arm vs agd exp arm
  res_AB$est <- summary(coxobj_dat_adj)$conf.int[1]
  mu <- summary(coxobj_dat_adj)$coef[1]
  sig <- summary(coxobj_dat_adj)$coef[4]
  res_AB$se <- sqrt((exp(sig^2) - 1) * exp(2 * mu + sig^2)) # log normal parametrization
  res_AB$ci_l <- summary(coxobj_dat_adj)$conf.int[3]
  res_AB$ci_u <- summary(coxobj_dat_adj)$conf.int[4]
  res_AB$pval <- summary(coxobj_dat_adj)$coef[6]

  # : get bootstrapped estimates if applicable
  if (!is.null(weights_object$boot)) {
    tmp_boot_obj <- weights_object$boot
    k <- dim(tmp_boot_obj)[3]

    cli::cli_progress_bar("Going through bootstrapped weights", total = k, .envir = .GlobalEnv)

    tmp_boot_est <- sapply(1:k, function(ii) {
      cli::cli_progress_update(.envir = .GlobalEnv)

      boot_x <- tmp_boot_obj[, , ii]
      boot_ipd_id <- weights_object$data$USUBJID[boot_x[, 1]]
      boot_ipd <- ipd[match(boot_ipd_id, ipd$USUBJID), , drop = FALSE]
      boot_ipd$weights <- boot_x[, 2]

      boot_dat <- rbind(boot_ipd, pseudo_ipd)
      boot_dat$ARM <- factor(boot_dat$ARM, levels = c(trt_agd, trt_ipd))

      boot_coxobj_dat_adj <- coxph(Surv(TIME, EVENT) ~ ARM, boot_dat, weights = weights)
      boot_AB_est <- summary(boot_coxobj_dat_adj)$coef[1]
      exp(boot_AB_est)
    })

    cli::cli_progress_done(.envir = .GlobalEnv)

    res$inferential[["boot_est"]] <- tmp_boot_est
  } else {
    res$inferential[["boot_est"]] <- NULL
  }

  # : make analysis report table
  res$inferential[["report_overall_robustCI"]] <- rbind(
    report_table_tte(coxobj_dat, medSurv_dat, tag = paste0("Before matching/", endpoint_name)),
    report_table_tte(coxobj_dat_adj, medSurv_dat_adj, tag = paste0("After matching/", endpoint_name))
  )

  if (is.null(res$inferential[["boot_est"]])) {
    res$inferential[["report_overall_bootCI"]] <- NULL
  } else {
    boot_res_AB <- res_AB
    boot_logres_se <- sd(log(res$inferential[["boot_est"]]), na.rm = TRUE)
    if (boot_ci_is_quantile) {
      boot_res_AB$ci_l <- quantile(res$inferential[["boot_est"]], p = 0.025)
      boot_res_AB$ci_u <- quantile(res$inferential[["boot_est"]], p = 0.975)
    } else {
      boot_res_AB$ci_l <- exp(log(boot_res_AB$est) + qnorm(0.025) * boot_logres_se)
      boot_res_AB$ci_u <- exp(log(boot_res_AB$est) + qnorm(0.975) * boot_logres_se)
    }
    tmp_report_table_tte <- report_table_tte(coxobj_dat_adj,
      medSurv_dat_adj,
      tag = paste0("After matching/", endpoint_name)
    )
    tmp_report_table_tte$`HR[95% CI]`[1] <- paste0(
      format(round(boot_res_AB$est, 2), nsmall = 2), "[",
      format(round(boot_res_AB$ci_l, 2), nsmall = 2), ";",
      format(round(boot_res_AB$ci_u, 2), nsmall = 2), "]"
    )
    res$inferential[["report_overall_bootCI"]] <- rbind(
      report_table_tte(coxobj_dat, medSurv_dat, tag = paste0("Before matching/", endpoint_name)),
      tmp_report_table_tte
    )
  }

  # output
  res
}

# MAIC inference functions for Binary outcome type ------------

maic_unanchored_binary <- function(res,
                                   res_AB,
                                   dat,
                                   ipd,
                                   pseudo_ipd,
                                   binary_robust_cov_type,
                                   weights_object,
                                   endpoint_name,
                                   eff_measure,
                                   boot_ci_is_quantile,
                                   trt_ipd,
                                   trt_agd) {
  # ~~~ Analysis table
  # : set up proper link
  glm_link <- switch(eff_measure,
    "RD" = poisson(link = "identity"),
    "RR" = poisson(link = "log"),
    "OR" = binomial(link = "logit")
  )

  # : fit glm for binary outcome and robust estimate with weights
  binobj_dat <- glm(RESPONSE ~ ARM, dat, family = glm_link)
  binobj_dat_adj <- glm(RESPONSE ~ ARM, dat, weights = weights, family = glm_link)
  bin_robust_cov <- clubSandwich::vcovCR(binobj_dat_adj,
    cluster = dat$USUBJID,
    type = binary_robust_cov_type
  )
  bin_robust_coef <- clubSandwich::conf_int(binobj_dat_adj, bin_robust_cov, coef = 2, p_values = TRUE)

  res$inferential[["model_before"]] <- binobj_dat
  res$inferential[["model_after"]] <- binobj_dat_adj

  mu <- bin_robust_coef$beta
  sig <- bin_robust_coef$SE
  res_AB$ci_l <- bin_robust_coef$CI_L
  res_AB$ci_u <- bin_robust_coef$CI_U
  res_AB$pval <- bin_robust_coef$p_val

  if (eff_measure %in% c("RR", "OR")) {
    res_AB$est <- exp(mu)
    res_AB$se <- sqrt((exp(sig^2) - 1) * exp(2 * mu + sig^2)) # log normal parameterization
    res_AB$ci_l <- exp(res_AB$ci_l)
    res_AB$ci_u <- exp(res_AB$ci_u)
  } else if (eff_measure == "RD") {
    res_AB$est <- mu * 100
    res_AB$se <- sig * 100
    res_AB$ci_l <- res_AB$ci_l * 100
    res_AB$ci_u <- res_AB$ci_u * 100
  }

  # : get bootstrapped estimates if applicable
  if (!is.null(weights_object$boot)) {
    tmp_boot_obj <- weights_object$boot
    k <- dim(tmp_boot_obj)[3]

    cli::cli_progress_bar("Going through bootstrapped weights", total = k, .envir = .GlobalEnv)

    tmp_boot_est <- sapply(1:k, function(ii) {
      cli::cli_progress_update(.envir = .GlobalEnv)

      boot_x <- tmp_boot_obj[, , ii]
      boot_ipd_id <- weights_object$data$USUBJID[boot_x[, 1]]
      boot_ipd <- ipd[match(boot_ipd_id, ipd$USUBJID), , drop = FALSE]
      boot_ipd$weights <- boot_x[, 2]

      boot_dat <- rbind(boot_ipd, pseudo_ipd)
      boot_dat$ARM <- factor(boot_dat$ARM, levels = c(trt_agd, trt_ipd))

      # does not matter use robust se or not, point estimate will not change and calculation would be faster
      boot_binobj_dat_adj <- glm(RESPONSE ~ ARM, boot_dat, weights = weights, family = glm_link)
      boot_bin_robust_cov <- clubSandwich::vcovCR(binobj_dat_adj,
        cluster = dat$USUBJID,
        type = binary_robust_cov_type
      )
      boot_bin_robust_coef <- clubSandwich::conf_int(boot_binobj_dat_adj, boot_bin_robust_cov, coef = 2)
      boot_AB_est <- boot_bin_robust_coef$beta
      if (eff_measure %in% c("RR", "OR")) {
        boot_AB_est <- exp(boot_AB_est)
      } else if (eff_measure == "RD") {
        boot_AB_est <- boot_AB_est * 100
      }
    })

    cli::cli_progress_done(.envir = .GlobalEnv)

    res$inferential[["boot_est"]] <- tmp_boot_est
  } else {
    res$inferential[["boot_est"]] <- NULL
  }

  # : make analysis report table
  res$inferential[["report_overall_robustCI"]] <- rbind(
    report_table_binary(
      binobj_dat,
      tag = paste0("Before matching/", endpoint_name),
      eff_measure = eff_measure
    ),
    report_table_binary(
      binobj_dat_adj,
      res_AB,
      tag = paste0("After matching/", endpoint_name),
      eff_measure = eff_measure
    )
  )

  if (is.null(res$inferential[["boot_est"]])) {
    res$inferential[["report_overall_bootCI"]] <- NULL
  } else {
    boot_res_AB <- res_AB
    boot_logres_se <- sd(log(res$inferential[["boot_est"]]), na.rm = TRUE)
    if (boot_ci_is_quantile) {
      boot_res_AB$ci_l <- quantile(res$inferential[["boot_est"]], p = 0.025)
      boot_res_AB$ci_u <- quantile(res$inferential[["boot_est"]], p = 0.975)
    } else {
      boot_res_AB$ci_l <- exp(log(boot_res_AB$est) + qnorm(0.025) * boot_logres_se)
      boot_res_AB$ci_u <- exp(log(boot_res_AB$est) + qnorm(0.975) * boot_logres_se)
    }
    tmp_report_table_binary <- report_table_binary(
      binobj_dat_adj,
      res_AB,
      tag = paste0("After matching/", endpoint_name),
      eff_measure = eff_measure
    )
    tmp_report_table_binary[[paste0(eff_measure, "[95% CI]")]][1] <- paste0(
      format(round(boot_res_AB$est, 2), nsmall = 2),
      "[",
      format(round(boot_res_AB$ci_l, 2), nsmall = 2),
      ";",
      format(round(boot_res_AB$ci_u, 2), nsmall = 2),
      "]"
    )
    res$inferential[["report_overall_bootCI"]] <- rbind(
      report_table_binary(
        binobj_dat,
        tag = paste0("Before matching/", endpoint_name),
        eff_measure = eff_measure
      ),
      tmp_report_table_binary
    )
  }
  # output
  res
}
