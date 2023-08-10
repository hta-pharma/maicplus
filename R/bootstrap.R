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
#' @param centered_colnames A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in internal dataset.
#' @param internal_time_name name of the time variable in ipd_centered (for time to event outcome)
#' @param internal_event_name name of the event variable in ipd_centered (for time to event outcome)
#' @param psuedo_ipd A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   Pseudo IPD is only used for unanchored case.
#' @param min_weight A numeric value that defines the minimum weight allowed.
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @examples
#' # Unanchored example
#' load(system.file("extdata", "ipd_centered.rda", package = "maicplus", mustWork = TRUE))
#' centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
#' centered_colnames <- paste0(centered_colnames, "_CENTERED")
#' pseudo_ipd <- read.csv(system.file("extdata", "psuedo_IPD.csv", package = "maicplus", mustWork = TRUE))
#' pseudo_ipd$ARM <- "B" # Need to specify ARM for pseudo ipd
#'
#' HR_bootstraps <- boot(
#'   data = ipd_centered,
#'   statistic = bootstrap_HR,
#'   centered_colnames = centered_colnames,
#'   internal_time_name = "TIME",
#'   internal_event_name = "EVENT",
#'   pseudo_ipd = pseudo_ipd,
#'   R = 1000
#' )
#'
#' @return The HR bootstraps
#' @export

bootstrap_HR <- function(ipd_centered, i, centered_colnames,
                         internal_time_name = "TIME", internal_event_name = "EVENT",
                         pseudo_ipd = NULL, min_weight = 0.0001) {
  # Resamples the centered internal IPD data
  bootstrap_data <- ipd_centered[i, ]

  # Estimates weights
  match_res <- estimate_weights(
    data = bootstrap_data,
    centered_colnames = centered_colnames
  )
  ipd_matched <- match_res$data

  if (length(unique(ipd_matched$ARM)) == 2) { # this is the anchored case
    combined_data <- ipd_matched
  } else if (length(unique(ipd_matched$ARM)) == 1) { # this is the unanchored case

    if (is.null(pseudo_ipd)) stop("For unanchored case pseudo_ipd needs to be specified")
    if (is.null(pseudo_ipd$ARM)) stop("external data should have an ARM column")

    combined_data <- merge_two_data(
      pseudo_ipd = pseudo_ipd,
      ipd_matched = ipd_matched,
      internal_time_name = "TIME",
      internal_event_name = "EVENT"
    )
  } else {
    stop("ipd_matched ARM should have 1 or 2 treatments")
  }

  # set weights that are below min_weight to min_weight to avoid issues with 0 values
  combined_data$weights <- ifelse(combined_data$weights < min_weight, min_weight, combined_data$weights)
  model <- as.formula(paste0("Surv(", internal_time_name, ", ", internal_event_name, "==1) ~ ARM"))
  cox_model <- survival::coxph(model, data = combined_data, weights = weights, robust = TRUE)
  HR <- exp(cox_model$coefficients)

  return(HR)
}
