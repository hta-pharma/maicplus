#' Bootstrapping weighted hazard ratios to obtain confidence interval
#'
#' This function uses \code{\link{boot}} to resample internal IPD data
#' to obtain bootstrapped hazard ratio values which is used to form confidence
#' interval. Note that due to lack of full IPD of covariates in the
#' aggregate data, pseudo IPD will not be bootstrapped.
#'
#' @param ipd_centered A data frame containing individual patient data
#'   from the internal IPD study. This data frame should have already been
#'   centered using [center_ipd] function
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param centered_colnames A character vector giving the names of the covariates to use in matching. These names must match the column names in internal dataset.
#' @param psuedo_ipd A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   Pseudo IPD is only used for unanchored case. pseudo_ipd should have time, event,
#'   and ARM variable with a same name and order as that of ipd_centered
#' @param model Survival model to run to calculate HR with correct time, event, and #'   ARM names
#' @param ref_treat Reference treatment for the survival model
#' @param min_weight A numeric value that defines the minimum weight allowed.
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @examples
#' # Unanchored example
#' load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
#' load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
#' ipd_centered <- center_ipd(ipd = adsl, agd = agd)
#'
#' pseudo_ipd <- read.csv(system.file("extdata", "psuedo_IPD.csv", package = "maicplus", mustWork = TRUE))
#' pseudo_ipd$ARM <- "B" # Need to specify ARM for pseudo ipd
#'
#' # Need to ensure pseudo_ipd and ipd_matched have same names
#' colnames(pseudo_ipd) <- c("TIME", "EVENT", "ARM")
#' \dontrun{
#' HR_bootstraps <- boot(
#'   data = ipd_centered,
#'   statistic = bootstrap_HR,
#'   centered_colnames = centered_colnames,
#'   pseudo_ipd = pseudo_ipd,
#'   model = Surv(TIME, EVENT == 1) ~ ARM,
#'   ref_treat = "B",
#'   R = 1000
#' )
#' }
#' @return The HR bootstraps
#' @export

bootstrap_HR <- function(ipd_centered, i, centered_colnames,
                         pseudo_ipd = NULL, model,
                         min_weight = 0.0001, ref_treat = NULL) {
  # Resamples the centered internal IPD data
  bootstrap_data <- ipd_centered[i, ]

  # Estimates weights
  weighted_data <- estimate_weights(
    data = bootstrap_data,
    centered_colnames = centered_colnames
  )

  ipd_matched <- weighted_data$data
  if (length(unique(ipd_matched$ARM)) == 2) { # this is the anchored case
    combined_data <- ipd_matched
  } else if (length(unique(ipd_matched$ARM)) == 1) { # this is the unanchored case
    if (is.null(pseudo_ipd)) stop("For unanchored case pseudo_ipd needs to be specified")
    if (is.null(pseudo_ipd$ARM)) stop("external data should have an ARM column")

    # Assign weight of 1 for pseudo_ipd if it not defined already
    pseudo_ipd$weights <- 1
    combined_data <- rbind(
      ipd_matched[, colnames(pseudo_ipd)],
      pseudo_ipd
    )
  } else {
    stop("ipd_matched ARM should have 1 or 2 treatments")
  }

  # set the base treatment
  if (is.null(ref_treat)) stop("Need to specify a value for ref_treat")
  combined_data$ARM <- stats::relevel(as.factor(combined_data$ARM), ref = ref_treat)

  # set weights that are below min_weight to min_weight to avoid issues with 0 values
  combined_data$weights <- ifelse(combined_data$weights < min_weight, min_weight, combined_data$weights)
  cox_model <- survival::coxph(model, data = combined_data, weights = weights, robust = TRUE)
  HR <- exp(cox_model$coefficients)

  return(HR)
}
