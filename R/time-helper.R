
# Create an environment for settings
settings_env <- new.env()

#' Get and Set Time Conversion Factors
#'
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
#' # Native time format is years
#' set_time_conversion(days = 1/365.25, weeks = 1/52.17857, months = 1/12, years = 1)
#'
#' # Native time format is days
#' set_time_conversion(days = 1, weeks = 7, months = 365.25/12, years = 365.25)
#'
set_time_conversion <- function(days = 1, weeks = 7, months = 365.25/12, years = 365.25) {
  settings_env$time_conversion <- c(days = days, weeks = weeks, months = months, years = years)
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
#' @param as A time scale to convert to
#'
#' @return Returns a numeric vector calculated from `times / get_time_conversion(factor = as)`
get_time_as <- function(times, as = c("days", "weeks", "months", "years")){
  if (!is.numeric(times)) stop('times arguments must be numeric')
  as <- match.arg(as)
  times / get_time_conversion(as)
}