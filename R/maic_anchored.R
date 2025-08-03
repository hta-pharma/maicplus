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
#' @param normalize_weights logical, default is \code{FALSE}. If \code{TRUE},
#'   \code{scaled_weights} (normalized weights) in \code{weights_object$data} will be used.
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
#' @example inst/examples/maic_anchored_ex.R
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
  if (normalize_weights) {
    ipd$weights <- weights_object$data$scaled_weights[match(ipd$USUBJID, weights_object$data$USUBJID)]
  } else {
    ipd$weights <- weights_object$data$weights[match(ipd$USUBJID, weights_object$data$USUBJID)]
  }
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
      dat,
      ipd,
      pseudo_ipd,
      km_conf_type,
      time_scale,
      weights_object,
      endpoint_name,
      normalize_weights,
      trt_ipd,
      trt_agd,
      boot_ci_type
    )
  } else if (endpoint_type == "binary") {
    maic_anchored_binary(
      res,
      res_BC = NULL,
      dat,
      ipd,
      pseudo_ipd,
      binary_robust_cov_type,
      weights_object,
      endpoint_name,
      normalize_weights,
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
                              dat,
                              ipd,
                              pseudo_ipd,
                              km_conf_type,
                              time_scale,
                              weights_object,
                              endpoint_name,
                              normalize_weights,
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
  medSurv_out <- cbind(medSurv_out[, 1:6],
    `events%` = medSurv_out$events * 100 / medSurv_out$n.max,
    medSurv_out[7:ncol(medSurv_out)]
  )
  medSurv_out <- cbind(trt_ind = c("C", "B", "A")[match(medSurv_out$treatment, levels(dat$ARM))], medSurv_out)

  res$descriptive[["summary"]] <- medSurv_out

  # ~~~ Analysis table (Cox model) before and after matching
  # fit PH Cox regression model
  coxobj_ipd <- coxph(Surv(TIME, EVENT) ~ ARM, ipd) # robust = TRUE or not makes a diff
  coxobj_ipd_adj <- coxph(Surv(TIME, EVENT) ~ ARM, ipd, weights = weights, robust = TRUE)
  coxobj_agd <- coxph(Surv(TIME, EVENT) ~ ARM, pseudo_ipd)

  # derive ipd exp arm vs agd exp arm via bucher
  res_AC_unadj <- as.list(summary(coxobj_ipd)$coef)[c(1, 3)] # est, se
  res_AC <- as.list(summary(coxobj_ipd_adj)$coef)[c(1, 4)] # est, robust se
  if (is.null(res_BC)) res_BC <- as.list(summary(coxobj_agd)$coef)[c(1, 3)] # est, se

  names(res_AC_unadj) <- names(res_AC) <- names(res_BC) <- c("est", "se")

  coxobj_ipd_summary <- summary(coxobj_ipd)
  res_AC_unadj$ci_l <- coxobj_ipd_summary$conf.int[3]
  res_AC_unadj$ci_u <- coxobj_ipd_summary$conf.int[4]
  res_AC_unadj$pval <- as.vector(coxobj_ipd_summary$waldtest[3])

  coxobj_ipd_adj_summary <- summary(coxobj_ipd_adj)
  res_AC$ci_l <- coxobj_ipd_adj_summary$conf.int[3]
  res_AC$ci_u <- coxobj_ipd_adj_summary$conf.int[4]
  res_AC$pval <- as.vector(coxobj_ipd_adj_summary$waldtest[3])

  coxobj_agd_summary <- summary(coxobj_agd)
  res_BC$ci_l <- coxobj_agd_summary$conf.int[3]
  res_BC$ci_u <- coxobj_agd_summary$conf.int[4]
  res_BC$pval <- as.vector(coxobj_agd_summary$waldtest[3])

  res_AB <- bucher(res_AC, res_BC, conf_lv = 0.95)
  res_AB_unadj <- bucher(res_AC_unadj, res_BC, conf_lv = 0.95)

  # : get bootstrapped estimates if applicable
  if (!is.null(weights_object$boot)) {
    keep_rows <- setdiff(seq_len(nrow(weights_object$data)), weights_object$rows_with_missing)
    boot_ipd_id <- weights_object$data[keep_rows, "USUBJID", drop = FALSE]

    boot_ipd <- merge(boot_ipd_id, ipd, by = "USUBJID", all.x = TRUE)
    if (nrow(boot_ipd) != nrow(boot_ipd_id)) stop("ipd has multiple observations for some patients")
    boot_ipd <- boot_ipd[match(boot_ipd$USUBJID, boot_ipd_id$USUBJID), ]

    stat_fun <- function(data, index, w_obj, normalize) {
      r <- dynGet("r", ifnotfound = NA) # Get bootstrap iteration
      if (!is.na(r)) {
        if (!all(index == w_obj$boot[, 1, r])) stop("Bootstrap and weight indices don't match")
        boot_ipd <- data[w_obj$boot[, 1, r], ]
        boot_ipd$weights <- w_obj$boot[, 2, r]

        if (normalize) {
          boot_ipd$weights <- ave(
            boot_ipd$weights,
            boot_ipd$ARM,
            FUN = function(w) w / mean(w, na.rm = TRUE)
          )
        }
      }
      boot_coxobj_dat_adj <- coxph(Surv(TIME, EVENT) ~ ARM, boot_ipd, weights = boot_ipd$weights, robust = TRUE)
      boot_res_AC <- list(est = coef(boot_coxobj_dat_adj)[1], se = sqrt(vcov(boot_coxobj_dat_adj)[1, 1]))
      # temp method to source in variance of BC in AgD via monte carlo, may be removed in future
      res_BC_mc <- res_BC
      res_BC_mc$est <- rnorm(1, mean = res_BC$est, sd = res_BC$se)
      boot_res_AB <- bucher(boot_res_AC, res_BC_mc, conf_lv = 0.95)
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
      normalize = normalize_weights,
      strata = weights_object$boot_strata
    )
    boot_ci <- boot.ci(boot_res, type = boot_ci_type, w_obj = weights_object, normalize = normalize_weights)

    l_u_index <- switch(boot_ci_type,
      "norm" = list(2, 3, "normal"),
      "basic" = list(4, 5, "basic"),
      "stud" = list(4, 5, "student"),
      "perc" = list(4, 5, "percent"),
      "bca" = list(4, 5, "bca")
    )

    # boot results for A v B, method 1 (maybe retired in future version)
    boot_res_AB <- list(
      est = as.vector(exp(boot_res$t0[1])),
      se = NA,
      ci_l = exp(boot_ci[[l_u_index[[3]]]][l_u_index[[1]]]),
      ci_u = exp(boot_ci[[l_u_index[[3]]]][l_u_index[[2]]]),
      pval = NA
    )

    # boot results for A v C
    boot_ci_ac <- boot.ci(boot_res, type = boot_ci_type, w_obj = weights_object, index = c(4, 6))
    boot_res_AC <- list(
      est = as.vector(exp(boot_res$t0[4])),
      se = NA,
      ci_l = exp(boot_ci_ac[[l_u_index[[3]]]][l_u_index[[1]]]),
      ci_u = exp(boot_ci_ac[[l_u_index[[3]]]][l_u_index[[2]]]),
      pval = NA
    )

    # boot results for A v B, method 2
    boot_res_AC2 <- list(
      est = as.vector(boot_res$t0[4]),
      se = NA,
      ci_l = boot_ci_ac[[l_u_index[[3]]]][l_u_index[[1]]],
      ci_u = boot_ci_ac[[l_u_index[[3]]]][l_u_index[[2]]],
      pval = NA
    )
    boot_res_AC2$se <- find_SE_from_CI(boot_res_AC2$ci_l, boot_res_AC2$ci_u, 0.95, log = FALSE)
    boot_res_AB2 <- bucher(boot_res_AC2, res_BC, conf_lv = 0.95)
    boot_res_AB2 <- list(
      est = exp(boot_res_AB2$est),
      se = NA,
      ci_l = exp(boot_res_AB2$ci_l),
      ci_u = exp(boot_res_AB2$ci_u),
      pval = NA
    )
  } else {
    boot_res <- NULL
    boot_res_AB <- NULL
    boot_res_AB2 <- NULL
    boot_res_AC <- NULL
  }

  # transform
  res_AB$est <- exp(res_AB$est)
  res_AB$ci_l <- exp(res_AB$ci_l)
  res_AB$ci_u <- exp(res_AB$ci_u)
  res_AB_unadj$est <- exp(res_AB_unadj$est)
  res_AB_unadj$ci_l <- exp(res_AB_unadj$ci_l)
  res_AB_unadj$ci_u <- exp(res_AB_unadj$ci_u)

  res_AC$est <- exp(res_AC$est)
  res_AC_unadj$est <- exp(res_AC_unadj$est)
  res_BC$est <- exp(res_BC$est)

  # : report all raw fitted obj
  res$inferential[["fit"]] <- list(
    km_before_ipd = kmobj_ipd,
    km_after_ipd = kmobj_ipd_adj,
    km_agd = kmobj_agd,
    model_before_ipd = coxobj_ipd,
    model_after_ipd = coxobj_ipd_adj,
    model_agd = coxobj_agd,
    res_AC = res_AC,
    res_AC_unadj = res_AC_unadj,
    res_BC = res_BC,
    res_AB = res_AB,
    res_AB_unadj = res_AB_unadj,
    boot_res = boot_res,
    boot_res_AC = boot_res_AC,
    boot_res_AB_mc = boot_res_AB,
    boot_res_AB = boot_res_AB2
  )

  # : compile HR result
  res$inferential[["summary"]] <- data.frame(
    case = c("AC", "adjusted_AC", "BC", "AB", "adjusted_AB"),
    HR = c(
      summary(coxobj_ipd)$conf.int[1],
      summary(coxobj_ipd_adj)$conf.int[1],
      summary(coxobj_agd)$conf.int[1],
      res_AB_unadj$est, res_AB$est
    ),
    LCL = c(
      summary(coxobj_ipd)$conf.int[3],
      summary(coxobj_ipd_adj)$conf.int[3],
      summary(coxobj_agd)$conf.int[3],
      res_AB_unadj$ci_l, res_AB$ci_l
    ),
    UCL = c(
      summary(coxobj_ipd)$conf.int[4],
      summary(coxobj_ipd_adj)$conf.int[4],
      summary(coxobj_agd)$conf.int[4],
      res_AB_unadj$ci_u, res_AB$ci_u
    ),
    pval = c(
      summary(coxobj_ipd)$waldtest[3],
      summary(coxobj_ipd_adj)$waldtest[3],
      summary(coxobj_agd)$waldtest[3],
      res_AB_unadj$pval, res_AB$pval
    )
  )

  # output
  res
}

# MAIC inference functions for Binary outcome type ------------
maic_anchored_binary <- function(res,
                                 res_BC = NULL,
                                 dat,
                                 ipd,
                                 pseudo_ipd,
                                 binary_robust_cov_type,
                                 weights_object,
                                 endpoint_name,
                                 normalize_weights,
                                 eff_measure,
                                 trt_ipd,
                                 trt_agd,
                                 boot_ci_type) {
  # ~~~ Analysis table
  # : set up proper link
  glm_link <- switch(eff_measure,
    "RD" = "identity",
    "RR" = "log",
    "OR" = "logit"
  )
  res_template <- list(
    est = NA,
    se = NA,
    ci_l = NA,
    ci_u = NA,
    pval = NA
  )

  # : fit glm for binary outcome and robust estimate with weights
  binobj_ipd <- glm(RESPONSE ~ ARM, ipd, family = binomial(link = glm_link))
  binobj_ipd_adj <- suppressWarnings(glm(RESPONSE ~ ARM, ipd, weights = weights, family = binomial(link = glm_link)))
  binobj_agd <- glm(RESPONSE ~ ARM, pseudo_ipd, family = binomial(link = glm_link))

  bin_robust_cov <- sandwich::vcovHC(binobj_ipd_adj, type = binary_robust_cov_type)
  bin_robust_coef <- lmtest::coeftest(binobj_ipd_adj, vcov. = bin_robust_cov)
  bin_robust_ci <- lmtest::coefci(binobj_ipd_adj, vcov. = bin_robust_cov)

  # : make general summary
  glmDesc_ipd <- glm_makeup(binobj_ipd, legend = "IPD, before matching", weighted = FALSE)
  glmDesc_ipd_adj <- glm_makeup(binobj_ipd_adj, legend = "IPD, after matching", weighted = TRUE)
  glmDesc_agd <- glm_makeup(binobj_agd, legend = "AgD, external", weighted = FALSE)
  glmDesc <- rbind(glmDesc_ipd, glmDesc_ipd_adj, glmDesc_agd)
  glmDesc <- cbind(trt_ind = c("C", "B", "A")[match(glmDesc$treatment, levels(dat$ARM))], glmDesc)
  rownames(glmDesc) <- NULL
  res$descriptive[["summary"]] <- glmDesc

  # derive ipd exp arm vs agd exp arm via bucher
  res_AC <- res_template
  res_AC$est <- bin_robust_coef[2, "Estimate"]
  res_AC$se <- bin_robust_coef[2, "Std. Error"]
  res_AC$ci_l <- bin_robust_ci[2, "2.5 %"]
  res_AC$ci_u <- bin_robust_ci[2, "97.5 %"]
  res_AC$pval <- bin_robust_coef[2, "Pr(>|z|)"]

  # unadjusted AC
  res_AC_unadj <- res_template
  res_AC_unadj$est <- summary(binobj_ipd)$coefficients[2, "Estimate"]
  res_AC_unadj$se <- summary(binobj_ipd)$coefficients[2, "Std. Error"]
  res_AC_unadj$ci_l <- confint.default(binobj_ipd)[2, "2.5 %"]
  res_AC_unadj$ci_u <- confint.default(binobj_ipd)[2, "97.5 %"]
  res_AC_unadj$pval <- summary(binobj_ipd)$coefficients[2, "Pr(>|z|)"]

  # BC
  if (is.null(res_BC)) {
    res_BC <- res_template
    res_BC$est <- summary(binobj_agd)$coefficients[2, "Estimate"]
    res_BC$se <- summary(binobj_agd)$coefficients[2, "Std. Error"]
    res_BC$ci_l <- confint.default(binobj_agd)[2, "2.5 %"]
    res_BC$ci_u <- confint.default(binobj_agd)[2, "97.5 %"]
    res_BC$pval <- summary(binobj_agd)$coefficients[2, "Pr(>|z|)"]
  }

  # derive AB
  res_AB <- bucher(res_AC, res_BC, conf_lv = 0.95)
  res_AB_unadj <- bucher(res_AC_unadj, res_BC, conf_lv = 0.95)

  # : get bootstrapped estimates if applicable
  if (!is.null(weights_object$boot)) {
    keep_rows <- setdiff(seq_len(nrow(weights_object$data)), weights_object$rows_with_missing)
    boot_ipd_id <- weights_object$data[keep_rows, "USUBJID", drop = FALSE]

    boot_ipd <- merge(boot_ipd_id, ipd, by = "USUBJID", all.x = TRUE)
    if (nrow(boot_ipd) != nrow(boot_ipd_id)) stop("ipd has multiple observations for some patients")
    boot_ipd <- boot_ipd[match(boot_ipd$USUBJID, boot_ipd_id$USUBJID), ]

    stat_fun <- function(data, index, w_obj, eff_measure, normalize) {
      r <- dynGet("r", ifnotfound = NA) # Get bootstrap iteration
      if (!is.na(r)) {
        if (!all(index == w_obj$boot[, 1, r])) stop("Bootstrap and weight indices don't match")
        boot_ipd <- data[w_obj$boot[, 1, r], ]
        boot_ipd$weights <- w_obj$boot[, 2, r]

        if (normalize) {
          boot_ipd$weights <- ave(
            boot_ipd$weights,
            boot_ipd$ARM,
            FUN = function(w) w / mean(w, na.rm = TRUE)
          )
        }
      }

      boot_binobj_dat_adj <- suppressWarnings(
        glm(RESPONSE ~ ARM, boot_ipd, weights = boot_ipd$weights, family = binomial(link = glm_link))
      )
      boot_AC_est <- coef(boot_binobj_dat_adj)[2]
      boot_AC_var <- vcov(boot_binobj_dat_adj)[2, 2]

      boot_res_AC <- list(est = boot_AC_est, se = sqrt(boot_AC_var))
      # temp method to source in variance of BC in AgD via monte carlo, may be removed in future
      res_BC_mc <- res_BC
      res_BC_mc$est <- rnorm(1, mean = res_BC$est, sd = res_BC$se)

      boot_res_AB <- bucher(boot_res_AC, res_BC_mc, conf_lv = 0.95)

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
      normalize = normalize_weights,
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

    # boot results for A v B, method 1 (maybe retired in future version)
    boot_res_AB <- list(
      est = as.vector(transform_estimate(boot_res$t0[1])),
      se = NA,
      ci_l = transform_estimate(boot_ci[[l_u_index[[3]]]][l_u_index[[1]]]),
      ci_u = transform_estimate(boot_ci[[l_u_index[[3]]]][l_u_index[[2]]]),
      pval = NA
    )

    # boot results for A v C
    boot_ci_ac <- boot.ci(boot_res, type = boot_ci_type, w_obj = weights_object, index = c(4, 6))
    boot_res_AC <- list(
      est = as.vector(transform_estimate(boot_res$t0[4])),
      se = NA,
      ci_l = transform_estimate(boot_ci_ac[[l_u_index[[3]]]][l_u_index[[1]]]),
      ci_u = transform_estimate(boot_ci_ac[[l_u_index[[3]]]][l_u_index[[2]]]),
      pval = NA
    )

    # boot results for A v B, method 2
    boot_res_AC2 <- list(
      est = as.vector(boot_res$t0[4]),
      se = NA,
      ci_l = boot_ci_ac[[l_u_index[[3]]]][l_u_index[[1]]],
      ci_u = boot_ci_ac[[l_u_index[[3]]]][l_u_index[[2]]],
      pval = NA
    )
    boot_res_AC2$se <- find_SE_from_CI(boot_res_AC2$ci_l, boot_res_AC2$ci_u, 0.95, log = FALSE)
    boot_res_AB2 <- bucher(boot_res_AC2, res_BC, conf_lv = 0.95)
    boot_res_AB2 <- list(
      est = transform_estimate(boot_res_AB2$est),
      se = NA,
      ci_l = transform_estimate(boot_res_AB2$ci_l),
      ci_u = transform_estimate(boot_res_AB2$ci_u),
      pval = NA
    )
  } else {
    boot_res_AC <- NULL
    boot_res_AB <- NULL
    boot_res_AB2 <- NULL
    boot_res <- NULL
  }

  # transform effect measures
  if (eff_measure %in% c("RR", "OR")) {
    res_AB <- transform_ratio(res_AB)
    res_AB_unadj <- transform_ratio(res_AB_unadj)
    res_AC <- transform_ratio(res_AC)
    res_AC_unadj <- transform_ratio(res_AC_unadj)
    res_BC <- transform_ratio(res_BC)
  } else if (eff_measure == "RD") {
    res_AB <- transform_absolute(res_AB)
    res_AB_unadj <- transform_absolute(res_AB_unadj)
    res_AC <- transform_absolute(res_AC)
    res_AC_unadj <- transform_absolute(res_AC_unadj)
    res_BC <- transform_absolute(res_BC)
  }


  # report all raw fitted obj
  res$inferential[["fit"]] <- list(
    model_before_ipd = binobj_ipd,
    model_after_ipd = binobj_ipd_adj,
    model_agd = binobj_agd,
    res_AC = res_AC,
    res_AC_unadj = res_AC_unadj,
    res_BC = res_BC,
    res_AB = res_AB,
    res_AB_unadj = res_AB_unadj,
    boot_res = boot_res,
    boot_res_AC = boot_res_AC,
    boot_res_AB_mc = boot_res_AB,
    boot_res_AB = boot_res_AB2
  )

  # compile binary effect estimate result
  res$inferential[["summary"]] <- data.frame(
    case = c("AC", "adjusted_AC", "BC", "AB", "adjusted_AB"),
    EST = c(
      res_AC_unadj$est,
      res_AC$est,
      res_BC$est,
      res_AB_unadj$est,
      res_AB$est
    ),
    LCL = c(
      res_AC_unadj$ci_l,
      res_AC$ci_l,
      res_BC$ci_l,
      res_AB_unadj$ci_l,
      res_AB$ci_l
    ),
    UCL = c(
      res_AC_unadj$ci_u,
      res_AC$ci_u,
      res_BC$ci_u,
      res_AB_unadj$ci_u,
      res_AB$ci_u
    ),
    pval = c(
      res_AC_unadj$pval,
      res_AC$pval,
      res_BC$pval,
      res_AB_unadj$pval,
      res_AB$pval
    )
  )
  names(res$inferential[["summary"]])[2] <- eff_measure

  # output
  res
}
