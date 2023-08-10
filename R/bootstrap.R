
#' Bootstrapping for MAIC weighted hazard ratios
#'
#' @param ipd_centered A data frame containing individual patient data
#'   from the internal IPD study. This data frame should have already been
#'   centered using [center_ipd] function
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param centered_colnames A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in internal dataset.
#' @param internal_time_name name of the time variable in ipd_centered (for time to event outcome)
#' @param internal_event_name how the event variable is named in the internal IPD (for time to event outcome)
#' @param external A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   The outcome variables names must match intervention_data.
#' @param min_weight A numeric value that defines the minimum weight allowed. 
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @details This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing MAIC weighting using
#'   {\link{estimate_weights}} and returning a weighted hazard ratio (HR) from a
#'   Cox proportional hazards model. This is used as the 'statistic' argument in
#'   the boot function.
#'
#' @return The HR as a numeric value.
#' @export

bootstrap_HR <- function(ipd_centered, i, centered_colnames, 
                         internal_time_name = "TIME", internal_event_name = "EVENT",
                         external, min_weight = 0.0001){
  
  # Resamples the centered internal IPD data
  bootstrap_data <- ipd_centered[i,]
  
  # Estimates weights
  match_res <- estimate_weights(data = bootstrap_data, 
                                centered_colnames = centered_colnames)

  internal_matched <- match_res$data
  if(is.null(external$ARM)){
    stop("external data should have an ARM column")
  }
  combined_data <- merge_two_data(external = external,
                                  internal = internal_matched, 
                                  internal_time_name = "TIME", 
                                  internal_event_name = "EVENT")
  
  if(length(unique(combined_data$ARM)) == 2){ #unanchored case

    # set weights that are below min_weight to min_weight to avoid issues with 0 values
    combined_data$weights <- ifelse(combined_data$weights < min_weight, min_weight, combined_data$weights)
    
    model <- as.formula(paste0("Surv(", internal_time_name, ", ", internal_event_name, "==1) ~ ARM"))
    cox_model <- survival::coxph(model, data = combined_data, weights = weights)
    HR <- exp(cox_model$coefficients)
  }
  return(HR)
}