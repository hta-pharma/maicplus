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
#' @param normalize_weights logical, default is \code{FALSE}. If \code{TRUE},
#'   \code{scaled_weights} (normalized weights) in \code{weights_object$data} will be used.
#' @param endpoint_type a string, one out of the following "binary", "tte" (time to event)
#' @param eff_measure a string, "RD" (risk difference), "OR" (odds ratio), "RR" (relative risk) for a binary endpoint;
#'   "HR" for a time-to-event endpoint. By default is \code{NULL}, "OR" is used for binary case, otherwise "HR" is used.
#' @param boot_ci_type a string, one of `c("norm","basic", "stud", "perc", "bca")` to select the type of bootstrap
#'   confidence interval. See [boot::boot.ci] for more details.
#' @param endpoint_name a string, name of time to event endpoint, to be show in the last line of title
#' @param time_scale a string, time unit of median survival time, taking a value of 'years', 'months', 'weeks' or
#'   'days'. NOTE: it is assumed that values in TIME column of \code{ipd} and \code{pseudo_ipd} is in the unit of days
#' @param km_conf_type a string, pass to \code{conf.type} of \code{survfit}
#' @param binary_robust_cov_type a string to pass to argument `type` of [sandwich::vcovHC], see possible options in the
#'   documentation of that function. Default is `"HC3"`
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
#' @importFrom survival survfit Surv coxph
#' @importFrom lmtest coeftest coefci
#' @importFrom sandwich vcovHC
#' @importFrom boot boot boot.ci
#' @return A list, contains 'descriptive' and 'inferential'
#' @example inst/examples/maic_unanchored_ex.R
#' @example inst/examples/maic_unanchored_binary_ex.R
#' @export

maic_unanchored <- function(weights_object,
                            ipd,
                            pseudo_ipd,
                            trt_ipd,
                            trt_agd,
                            trt_var_ipd = "ARM",
                            trt_var_agd = "ARM",
                            normalize_weights = FALSE,
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

  res_AB_unadj <- res_AB <- list(
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
    binary_robust_cov_type <- match.arg(
      binary_robust_cov_type,
      choices = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5")
    )
  } else if (endpoint_type == "tte") { # for time to event effect measure

    if (!all(c("USUBJID", "TIME", "EVENT", trt_var_ipd) %in% names(ipd))) {
      stop("ipd needs to include at least USUBJID, TIME, EVENT, ", trt_var_ipd)
    }
    if (!all(c("TIME", "EVENT", trt_var_agd) %in% names(pseudo_ipd))) {
      stop("pseudo_ipd needs to include at least TIME, EVENT, ", trt_var_agd)
    }
    eff_measure <- match.arg(eff_measure, choices = c("HR"), several.ok = FALSE)
  }
  boot_ci_type <- match.arg(boot_ci_type)

  # ==> IPD and AgD data preparation ------------------------------------------
  # : subset ipd, retain only ipd from interested trts
  ipd <- ipd[ipd$ARM == trt_ipd, , drop = TRUE]
  pseudo_ipd <- pseudo_ipd[pseudo_ipd$ARM == trt_agd, , drop = TRUE]

  # : assign weights to real and pseudo ipd
  if (normalize_weights) {
    ipd$weights <- weights_object$data$scaled_weights[match(weights_object$data$USUBJID, ipd$USUBJID)]
  } else {
    ipd$weights <- weights_object$data$weights[match(weights_object$data$USUBJID, ipd$USUBJID)]
  }
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
      res, res_AB, res_AB_unadj, dat, ipd, pseudo_ipd, km_conf_type, time_scale,
      weights_object, endpoint_name, normalize_weights, boot_ci_type, trt_ipd, trt_agd
    )
  } else if (endpoint_type == "binary") {
    maic_unanchored_binary(
      res, res_AB, res_AB_unadj, dat, ipd, pseudo_ipd, binary_robust_cov_type,
      weights_object, endpoint_name, normalize_weights, eff_measure, boot_ci_type, trt_ipd, trt_agd
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
                                res_AB_unadj,
                                dat,
                                ipd,
                                pseudo_ipd,
                                km_conf_type,
                                time_scale,
                                weights_object,
                                endpoint_name,
                                normalize_weights,
                                boot_ci_type,
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
  medSurv_out <- cbind(trt_ind = c("B", "A")[match(medSurv_out$treatment, levels(dat$ARM))], medSurv_out)

  res$descriptive[["summary"]] <- medSurv_out

  # ~~~ Analysis table (Cox model) before and after matching
  # : fit PH Cox regression model
  coxobj_dat <- coxph(Surv(TIME, EVENT) ~ ARM, dat)
  coxobj_dat_adj <- coxph(Surv(TIME, EVENT) ~ ARM, dat, weights = weights, robust = TRUE)

  # : derive adjusted estimate for ipd exp arm vs agd exp arm
  res_AB$est <- summary(coxobj_dat_adj)$conf.int[1]
  mu <- summary(coxobj_dat_adj)$coef[1]
  sig <- summary(coxobj_dat_adj)$coef[4]
  res_AB$se <- sqrt((exp(sig^2) - 1) * exp(2 * mu + sig^2)) # log normal parametrization
  res_AB$ci_l <- summary(coxobj_dat_adj)$conf.int[3]
  res_AB$ci_u <- summary(coxobj_dat_adj)$conf.int[4]
  res_AB$pval <- summary(coxobj_dat_adj)$coef[6]

  # : derive unadjusted estimate
  res_AB_unadj$est <- summary(coxobj_dat)$conf.int[1]
  mu <- summary(coxobj_dat)$coef[1]
  sig <- summary(coxobj_dat)$coef[3]
  res_AB_unadj$se <- sqrt((exp(sig^2) - 1) * exp(2 * mu + sig^2)) # log normal parametrization
  res_AB_unadj$ci_l <- summary(coxobj_dat)$conf.int[3]
  res_AB_unadj$ci_u <- summary(coxobj_dat)$conf.int[4]
  res_AB_unadj$pval <- summary(coxobj_dat)$coef[5]

  # : get bootstrapped estimates if applicable
  if (!is.null(weights_object$boot)) {
    keep_rows <- setdiff(seq_len(nrow(weights_object$data)), weights_object$rows_with_missing)
    boot_ipd_id <- weights_object$data[keep_rows, "USUBJID", drop = FALSE]

    boot_ipd <- merge(boot_ipd_id, ipd, by = "USUBJID", all.x = TRUE)
    if (nrow(boot_ipd) != nrow(boot_ipd_id)) stop("ipd has multiple observations for some patients")
    boot_ipd <- boot_ipd[match(boot_ipd$USUBJID, boot_ipd_id$USUBJID), ]

    stat_fun <- function(data, index, w_obj, pseudo_ipd, normalize) {
      boot_ipd <- data[index, ]
      r <- dynGet("r", ifnotfound = NA) # Get bootstrap iteration
      if (!is.na(r)) {
        if (!all(index == w_obj$boot[, 1, r])) stop("Bootstrap and weight indices don't match")
        boot_ipd$weights <- w_obj$boot[, 2, r]
        if (normalize) boot_ipd$weights <- boot_ipd$weights / mean(boot_ipd$weights, na.rm = TRUE)
      }
      boot_dat <- rbind(boot_ipd, pseudo_ipd)
      boot_dat$ARM <- factor(boot_dat$ARM, levels = c(trt_agd, trt_ipd))
      boot_coxobj_dat_adj <- coxph(Surv(TIME, EVENT) ~ ARM, boot_dat, weights = weights)
      c(est = coef(boot_coxobj_dat_adj)[1], var = vcov(boot_coxobj_dat_adj)[1, 1])
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
      pseudo_ipd = pseudo_ipd,
      normalize = normalize_weights,
      strata = weights_object$boot_strata
    )
    boot_ci <- boot.ci(boot_res, type = boot_ci_type, w_obj = weights_object, pseudo_ipd = pseudo_ipd)

    l_u_index <- switch(boot_ci_type,
      "norm" = list(2, 3, "normal"),
      "basic" = list(4, 5, "basic"),
      "stud" = list(4, 5, "student"),
      "perc" = list(4, 5, "percent"),
      "bca" = list(4, 5, "bca")
    )

    transform_estimate <- exp
    boot_res_AB <- list(
      est = as.vector(transform_estimate(boot_res$t0[1])),
      se = NA,
      ci_l = transform_estimate(boot_ci[[l_u_index[[3]]]][l_u_index[[1]]]),
      ci_u = transform_estimate(boot_ci[[l_u_index[[3]]]][l_u_index[[2]]]),
      pval = NA
    )
  } else {
    boot_res_AB <- NULL
    boot_res <- NULL
  }

  # : report all raw fitted obj
  res$inferential[["fit"]] <- list(
    km_before = kmobj_dat,
    km_after = kmobj_dat_adj,
    model_before = coxobj_dat,
    model_after = coxobj_dat_adj,
    res_AB = res_AB,
    res_AB_unadj = res_AB_unadj,
    boot_res = boot_res,
    boot_res_AB = boot_res_AB
  )

  # : compile HR result
  res$inferential[["summary"]] <- data.frame(
    case = c("AB", "adjusted_AB"),
    HR = c(res_AB_unadj$est, res_AB$est),
    LCL = c(res_AB_unadj$ci_l, res_AB$ci_l),
    UCL = c(res_AB_unadj$ci_u, res_AB$ci_u),
    pval = c(res_AB_unadj$pval, res_AB$pval)
  )

  # output
  res
}

# MAIC inference functions for Binary outcome type ------------

maic_unanchored_binary <- function(res,
                                   res_AB,
                                   res_AB_unadj,
                                   dat,
                                   ipd,
                                   pseudo_ipd,
                                   binary_robust_cov_type,
                                   weights_object,
                                   endpoint_name,
                                   normalize_weights,
                                   eff_measure,
                                   boot_ci_type,
                                   trt_ipd,
                                   trt_agd) {
  # ~~~ Analysis table
  # : set up proper link
  glm_link <- switch(eff_measure,
    "RD" = poisson(link = "identity"),
    "RR" = poisson(link = "log"),
    "OR" = binomial(link = "logit")
  )
  transform_estimate <- switch(eff_measure,
    "RD" = function(x) x * 100,
    "RR" = exp,
    "OR" = exp
  )

  # : fit glm for binary outcome and robust estimate with weights
  binobj_dat <- glm(RESPONSE ~ ARM, dat, family = glm_link)
  binobj_dat_adj <- glm(RESPONSE ~ ARM, dat, weights = weights, family = glm_link) |> suppressWarnings()

  bin_robust_cov <- sandwich::vcovHC(binobj_dat_adj, type = binary_robust_cov_type)
  bin_robust_coef <- lmtest::coeftest(binobj_dat_adj, vcov. = bin_robust_cov)
  bin_robust_ci <- lmtest::coefci(binobj_dat_adj, vcov. = bin_robust_cov)

  # : make general summary
  glmDesc_dat <- glm_makeup(binobj_dat, legend = "Before matching", weighted = FALSE)
  glmDesc_dat_adj <- glm_makeup(binobj_dat_adj, legend = "After matching", weighted = TRUE)
  glmDesc <- rbind(glmDesc_dat, glmDesc_dat_adj)
  glmDesc <- cbind(trt_ind = c("B", "A")[match(glmDesc$treatment, levels(dat$ARM))], glmDesc)
  rownames(glmDesc) <- NULL
  res$descriptive[["summary"]] <- glmDesc

  # : derive adjusted estimate
  res_AB$est <- bin_robust_coef[2, "Estimate"]
  res_AB$se <- bin_robust_coef[2, "Std. Error"]
  res_AB$ci_l <- bin_robust_ci[2, "2.5 %"]
  res_AB$ci_u <- bin_robust_ci[2, "97.5 %"]
  res_AB$pval <- bin_robust_coef[2, "Pr(>|z|)"]

  # : derive unadjusted estimate
  binobj_dat_summary <- summary(binobj_dat)
  res_AB_unadj$est <- binobj_dat_summary$coefficients[2, "Estimate"]
  res_AB_unadj$se <- binobj_dat_summary$coefficients[2, "Std. Error"]
  res_AB_unadj$ci_l <- confint.default(binobj_dat)[2, "2.5 %"]
  res_AB_unadj$ci_u <- confint.default(binobj_dat)[2, "97.5 %"]
  res_AB_unadj$pval <- binobj_dat_summary$coefficients[2, "Pr(>|z|)"]

  # : transform
  if (eff_measure %in% c("RR", "OR")) {
    res_AB <- transform_ratio(res_AB)
    res_AB_unadj <- transform_ratio(res_AB_unadj)
  } else if (eff_measure == "RD") {
    res_AB <- transform_absolute(res_AB)
    res_AB_unadj <- transform_absolute(res_AB_unadj)
  }

  # : get bootstrapped estimates if applicable
  if (!is.null(weights_object$boot)) {
    keep_rows <- setdiff(seq_len(nrow(weights_object$data)), weights_object$rows_with_missing)
    boot_ipd_id <- weights_object$data[keep_rows, "USUBJID", drop = FALSE]

    boot_ipd <- merge(boot_ipd_id, ipd, by = "USUBJID", all.x = TRUE)
    if (nrow(boot_ipd) != nrow(boot_ipd_id)) stop("ipd has multiple observations for some patients")
    boot_ipd <- boot_ipd[match(boot_ipd$USUBJID, boot_ipd_id$USUBJID), ]

    stat_fun <- function(data, index, w_obj, pseudo_ipd, normalize) {
      boot_ipd <- data[index, ]
      r <- dynGet("r", ifnotfound = NA) # Get bootstrap iteration
      if (!is.na(r)) {
        if (!all(index == w_obj$boot[, 1, r])) stop("Bootstrap and weight indices don't match")
        boot_ipd$weights <- w_obj$boot[, 2, r]
        if (normalize) boot_ipd$weights <- boot_ipd$weights / mean(boot_ipd$weights, na.rm = TRUE)
      }
      boot_dat <- rbind(boot_ipd, pseudo_ipd)
      boot_dat$ARM <- factor(boot_dat$ARM, levels = c(trt_agd, trt_ipd))
      boot_binobj_dat_adj <- glm(RESPONSE ~ ARM, boot_dat, weights = weights, family = glm_link) |> suppressWarnings()
      c(est = coef(boot_binobj_dat_adj)[2], var = vcov(boot_binobj_dat_adj)[2, 2])
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
      pseudo_ipd = pseudo_ipd,
      normalize = normalize_weights,
      strata = weights_object$boot_strata
    )
    boot_ci <- boot.ci(boot_res, type = boot_ci_type, w_obj = weights_object, pseudo_ipd = pseudo_ipd)

    l_u_index <- switch(boot_ci_type,
      "norm" = list(2, 3, "normal"),
      "basic" = list(4, 5, "basic"),
      "stud" = list(4, 5, "student"),
      "perc" = list(4, 5, "percent"),
      "bca" = list(4, 5, "bca")
    )

    boot_res_AB <- list(
      est = as.vector(transform_estimate(boot_res$t0[1])),
      se = NA,
      ci_l = transform_estimate(boot_ci[[l_u_index[[3]]]][l_u_index[[1]]]),
      ci_u = transform_estimate(boot_ci[[l_u_index[[3]]]][l_u_index[[2]]]),
      pval = NA
    )
  } else {
    boot_res_AB <- NULL
    boot_res <- NULL
  }

  # : report all raw fitted obj
  res$inferential[["fit"]] <- list(
    model_before = binobj_dat,
    model_after = binobj_dat_adj,
    res_AB = res_AB,
    res_AB_unadj = res_AB_unadj,
    boot_res = boot_res,
    boot_res_AB = boot_res_AB
  )

  # : compile binary effect estimate result
  res$inferential[["summary"]] <- data.frame(
    case = c("AB", "adjusted_AB"),
    EST = c(
      res_AB_unadj$est,
      res_AB$est
    ),
    LCL = c(
      res_AB_unadj$ci_l,
      res_AB$ci_l
    ),
    UCL = c(
      res_AB_unadj$ci_u,
      res_AB$ci_u
    ),
    pval = c(
      res_AB_unadj$pval,
      res_AB$pval
    )
  )
  names(res$inferential[["summary"]])[2] <- eff_measure

  # : output
  res
}
