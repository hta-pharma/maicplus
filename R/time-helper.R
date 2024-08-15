# Create an environment for settings
settings_env <- new.env()

#' Get and Set Time Conversion Factors
#'
#' @param default The default time scale, commonly whichever has factor = 1
#' @param days Factor to divide data time units to get time in days
#' @param weeks Factor to divide data time units to get time in weeks
#' @param months Factor to divide data time units to get time in months
#' @param years Factor to divide data time units to get time in years
#'
#' @return No value returned. Conversion factors are stored internally and used within functions.
#' @export
#' @rdname time_conversion
#'
#' @examples
#' # The default time scale is days:
#' set_time_conversion(default = "days", days = 1, weeks = 7, months = 365.25 / 12, years = 365.25)
#'
#' # Set the default time scale to years
#' set_time_conversion(
#'   default = "years",
#'   days = 1 / 365.25,
#'   weeks = 1 / 52.17857,
#'   months = 1 / 12,
#'   years = 1
#' )
#'
set_time_conversion <- function(default = "days", days = 1, weeks = 7, months = 365.25 / 12, years = 365.25) {
  if (!default %in% c("days", "weeks", "months", "years")) {
    stop("default must be one of \"days\", \"weeks\", \"months\", \"years\")")
  }
  factors <- c(days = days, weeks = weeks, months = months, years = years)
  check_factors <- vapply(factors, function(x) isFALSE(!is.finite(x) || x == 0), logical(1L))
  if (!all(check_factors)) {
    stop(
      "Conversion factors must be finite non-zero numerical values: ",
      paste0(names(factors)[!check_factors], " = ", factors[!check_factors], collapse = ", ")
    )
  }
  settings_env$time_conversion <- factors
  settings_env$default_time_scale <- default
}


#' @param factor Time factor to get.
#' @rdname time_conversion
#' @export
#'
#' @examples
#' # Get time scale factors:
#' get_time_conversion("years")
#' get_time_conversion("weeks")
get_time_conversion <- function(factor = c("days", "weeks", "months", "years")) {
  factor <- match.arg(factor, several.ok = TRUE)
  if (!exists("time_conversion", settings_env)) {
    warning("No time conversion factors previously set. Setting defaults.")
    set_time_conversion()
  }
  settings_env$time_conversion[factor]
}


#' Convert Time Values Using Scaling Factors
#'
#' @param times Numeric time values
#' @param as A time scale to convert to. One of "days", "weeks", "months", "years"
#'
#' @return Returns a numeric vector calculated from `times / get_time_conversion(factor = as)`
#' @export
#' @examples
#' get_time_as(50, as = "months")
get_time_as <- function(times, as = NULL) {
  if (is.null(as)) as <- settings_env$default_time_scale
  if (!is.numeric(times)) stop("times arguments must be numeric")
  as <- match.arg(as, c("days", "weeks", "months", "years"))
  times / get_time_conversion(as)
}
