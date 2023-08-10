#' Bootstrapping for MAIC weighted hazard ratios
#'
#' @param intervention_data A data frame containing individual patient data from the internal IPD study including columns:
#' treatment, time, status, centered baseline characteristics.
#' @param centered_colnames A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param model A model formula in the form 'Surv(Time, Event==1) ~ ARM'.
#'   Variable names need to match the corresponding columns in intervention_data.
#' @param comparator_data A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect
#'   including columns: time, event, weight=1, scaled_weight=1.
#'   The outcome variables names must match intervention_data.
#' @param min_weight A numeric value that defines the minimum weight allowed.
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#' @param trt_ext A character giving the name of the comparator treatment to be used as reference
#'
#' @details This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing MAIC weighting using
#'   {\link{estimate_weights}} and returning a weighted hazard ratio (HR) from a
#'   Cox proportional hazards model. This is used as the 'statistic' argument in
#'   the boot function.
#'
#' @return The HR as a numeric value.
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#' @example inst/examples/MAIC_example_analysis.R
#' @export



bootstrap_HR <- function(intervention_data, centered_colnames, i, model, comparator_data,  min_weight = 0.0001, trt_ext){

  # create a visible binding for R CMD check
  wt <- NULL

  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(data=bootstrap_data,  centered_colnames=centered_colnames)


  # Give comparator data weights of 1
  comparator_data_wts <- comparator_data %>% dplyr::mutate(weights=1, scaled_weights=1)

  # Add the comparator data
  combined_data <- dplyr::bind_rows(perform_wt$data, comparator_data_wts)
  combined_data$treatment <- stats::relevel(as.factor(combined_data$treatment),ref=trt_ext)

  # set weights that are below min_weight to min_weight to avoid issues with 0 values
  combined_data$wt <- ifelse(combined_data$weights < min_weight, min_weight, combined_data$weights)

  # survival data stat
  cox_model <- survival::coxph(model, data = combined_data, weights = wt)
  HR <- exp(cox_model$coefficients)
}



