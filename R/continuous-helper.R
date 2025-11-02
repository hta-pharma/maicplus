#' Create pseudo IPD given aggregated continuous data
#'
#' @param continuous_agd a data.frame with columns ARM, MEAN, SD, N representing aggregated
#'   continuous outcome data
#'
#' @return a data.frame of pseudo continuous IPD, with columns USUBJID, ARM, AVAL
#' @example inst/examples/get_pseudo_ipd_continuous_ex.R
#' @export

get_pseudo_ipd_continuous <- function(continuous_agd) {
  # pre-check
  if (!is.data.frame(continuous_agd)) {
    stop("continuous_agd should be a data.frame with columns 'ARM', 'MEAN', 'SD', 'N'")
  }

  names(continuous_agd) <- toupper(names(continuous_agd))

  if (!all(c("ARM", "MEAN", "SD", "N") %in% names(continuous_agd))) {
    stop("continuous_agd should be a data.frame with columns 'ARM', 'MEAN', 'SD', 'N'")
  }

  if (!is.numeric(continuous_agd$MEAN) || !is.numeric(continuous_agd$SD) || !is.numeric(continuous_agd$N)) {
    stop("MEAN, SD, and N columns must be numeric")
  }

  if (any(continuous_agd$SD < 0)) {
    stop("SD values must be non-negative")
  }

  if (any(continuous_agd$N <= 0) || any(continuous_agd$N != floor(continuous_agd$N))) {
    stop("N values must be positive integers")
  }

  # create pseudo continuous IPD
  continuous_agd$ARM <- factor(continuous_agd$ARM, levels = unique(continuous_agd$ARM))

  # generate pseudo IPD for each arm
  pseudo_ipd_list <- lapply(seq_len(nrow(continuous_agd)), function(i) {
    arm <- continuous_agd$ARM[i]
    n <- continuous_agd$N[i]
    mean_val <- continuous_agd$MEAN[i]
    sd_val <- continuous_agd$SD[i]

    # generate values assuming normal distribution
    # for reproducibility in examples, we generate values that preserve mean and SD
    if (sd_val == 0) {
      # all values are the same
      aval <- rep(mean_val, n)
    } else {
      # generate standard normal, then scale to match mean and SD
      set.seed(as.integer(charToRaw(paste0(arm, i)))[1]) # deterministic per arm
      aval <- rnorm(n, mean = mean_val, sd = sd_val)
    }

    data.frame(
      USUBJID = paste0("pseudo_cont_subj_", arm, "_", seq_len(n)),
      ARM = as.character(arm),
      AVAL = aval,
      stringsAsFactors = FALSE
    )
  })

  # combine all arms
  pseudo_ipd <- do.call(rbind, pseudo_ipd_list)
  rownames(pseudo_ipd) <- NULL

  # output
  pseudo_ipd
}


#' Helper function to summarize outputs from lm fit
#'
#' @param lmobj returned object from \code{stats::lm}
#' @param legend label to indicate the linear model fit
#' @param weighted logical flag indicating whether weights have been applied in the lm fit
#' @returns A `data.frame` containing a summary of the number of subjects and descriptive
#'   statistics in a linear regression model.
#' @examples
#' data(adlb_sat)
#' pseudo_adlb <- get_pseudo_ipd_continuous(
#'   continuous_agd = data.frame(
#'     ARM = "B",
#'     MEAN = 50,
#'     SD = 10,
#'     N = 100
#'   )
#' )
#' combined_data <- rbind(adlb_sat[, c("USUBJID", "ARM", "AVAL")], pseudo_adlb)
#' combined_data$ARM <- as.factor(combined_data$ARM)
#' lmobj_dat <- stats::lm(AVAL ~ ARM, combined_data)
#' lm_makeup(lmobj_dat)
#' @export
lm_makeup <- function(lmobj, legend = "before matching", weighted = FALSE) {
  arm <- levels(lmobj$model$ARM)
  if (!weighted) {
    n <- tapply(seq_len(nrow(lmobj$model)), lmobj$model$ARM, length)
    mean_val <- tapply(lmobj$model$AVAL, lmobj$model$ARM, mean)
    sd_val <- tapply(lmobj$model$AVAL, lmobj$model$ARM, sd)
  } else {
    # get weights from model
    wts <- lmobj$model$`(weights)`
    n <- tapply(wts, lmobj$model$ARM, sum)
    # weighted mean and SD
    mean_val <- tapply(seq_along(lmobj$model$AVAL), lmobj$model$ARM, function(idx) {
      stats::weighted.mean(lmobj$model$AVAL[idx], wts[idx])
    })
    sd_val <- tapply(seq_along(lmobj$model$AVAL), lmobj$model$ARM, function(idx) {
      aval <- lmobj$model$AVAL[idx]
      w <- wts[idx]
      weighted_mean <- stats::weighted.mean(aval, w)
      sqrt(sum(w * (aval - weighted_mean)^2) / sum(w))
    })
  }

  data.frame(
    treatment = arm,
    type = legend,
    n = n,
    mean = mean_val,
    sd = sd_val
  )
}


#' Calculate standardized mean difference (SMD)
#'
#' @param mean1 Mean of group 1
#' @param mean2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param method Character string specifying SMD calculation method: "cohen" (Cohen's d),
#'   "hedges" (Hedges' g), or "glass" (Glass's delta)
#'
#' @return A numeric value representing the SMD
#' @keywords internal
calculate_smd <- function(mean1, mean2, sd1, sd2, n1, n2, method = "cohen") {
  method <- match.arg(method, c("cohen", "hedges", "glass"))

  mean_diff <- mean1 - mean2

  smd <- switch(method,
    "cohen" = {
      # Cohen's d: uses pooled standard deviation
      pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
      mean_diff / pooled_sd
    },
    "hedges" = {
      # Hedges' g: bias-corrected Cohen's d
      pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
      cohens_d <- mean_diff / pooled_sd
      correction_factor <- 1 - (3 / (4 * (n1 + n2 - 2) - 1))
      cohens_d * correction_factor
    },
    "glass" = {
      # Glass's delta: uses control group (group 2) SD only
      mean_diff / sd2
    }
  )

  return(smd)
}


#' Calculate standard error of standardized mean difference
#'
#' @param smd The standardized mean difference
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param method Character string specifying SMD calculation method
#'
#' @return A numeric value representing the SE of SMD
#' @keywords internal
calculate_smd_se <- function(smd, n1, n2, method = "cohen") {
  method <- match.arg(method, c("cohen", "hedges", "glass"))

  se <- switch(method,
    "cohen" = {
      # SE for Cohen's d
      sqrt((n1 + n2) / (n1 * n2) + smd^2 / (2 * (n1 + n2)))
    },
    "hedges" = {
      # SE for Hedges' g (approximately same as Cohen's d)
      correction_factor <- 1 - (3 / (4 * (n1 + n2 - 2) - 1))
      sqrt((n1 + n2) / (n1 * n2) + smd^2 / (2 * (n1 + n2))) * correction_factor
    },
    "glass" = {
      # SE for Glass's delta
      sqrt(n1 / (n1 * n2) + smd^2 / (2 * n2))
    }
  )

  return(se)
}

