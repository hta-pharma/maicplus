# Functions for matching step: estimation of individual weights

# functions to be exported ---------------------------------------

#' Derive individual weights in the matching step of MAIC
#'
#' Assuming data is properly processed, this function takes individual patient data (IPD) with centered covariates
#' (effect modifiers and/or prognostic variables) as input, and generates weights for each individual in IPD trial that
#' matches the chosen statistics of those covariates in Aggregated Data (AgD) trial.
#'
#' @param data a numeric matrix, centered covariates of IPD, no missing value in any cell is allowed
#' @param centered_colnames a character or numeric vector (column indicators) of centered covariates
#' @param start_val a scalar, the starting value for all coefficients of the propensity score regression
#' @param method a string, name of the optimization algorithm (see 'method' argument of \code{base::optim()}).
#' The default is `"BFGS"`, other options are `"Nelder-Mead"`, `"CG"`, `"L-BFGS-B"`, `"SANN"`, and `"Brent"`
#' @param ... all other arguments from \code{base::optim()}
#'
#' @return a list with the following 4 elements,
#' \describe{
#'   \item{data}{a data.frame, includes the input \code{data} with appended column 'weights' and 'scaled_weights'.
#'   Scaled weights has a summation to be the number of rows in \code{data} that has no missing value in any of the
#'   effect modifiers}
#'   \item{centered_colnames}{column names of centered effect modifiers in \code{data}}
#'   \item{nr_missing}{number of rows in \code{data} that has at least 1 missing value in specified centered effect
#'   modifiers}
#'   \item{ess}{effective sample size, square of sum divided by sum of squares}
#'   \item{opt}{R object returned by \code{base::optim()}, for assess convergence and other details}
#' }
#' @export

estimate_weights <- function(data, centered_colnames = NULL, start_val = 0, method = "BFGS", ...) {
  # pre check
  ch1 <- is.data.frame(data)
  if (!ch1) stop("'data' is not a data.frame")

  ch2 <- (!is.null(centered_colnames))
  if (ch2 && is.numeric(centered_colnames)) {
    ch2b <- any(centered_colnames < 1 | centered_colnames > ncol(data))
    if (ch2b) stop("specified centered_colnames are out of bound")
  } else if (ch2 && is.character(centered_colnames)) {
    ch2b <- !all(centered_colnames %in% names(data))
    if (ch2b) stop("1 or more specified centered_colnames are not found in 'data'")
  } else {
    stop("'centered_colnames' should be either a numeric or character vector")
  }

  ch3 <- sapply(seq_along(centered_colnames), function(ii) {
    !is.numeric(data[, centered_colnames[ii]])
  })
  if (any(ch3)) {
    stop(paste0("following columns of 'data' are not numeric for the calculation:", paste(which(ch3), collapse = ",")))
  }

  # prepare data for optimization
  if (is.null(centered_colnames)) centered_colnames <- seq_len(ncol(data))
  EM <- data[, centered_colnames, drop = FALSE]
  ind <- apply(EM, 1, function(xx) any(is.na(xx)))
  nr_missing <- sum(ind)
  EM <- as.matrix(EM[!ind, , drop = FALSE])

  # objective and gradient functions
  objfn <- function(alpha, X) {
    sum(exp(X %*% alpha))
  }
  gradfn <- function(alpha, X) {
    colSums(sweep(X, 1, exp(X %*% alpha), "*"))
  }

  # estimate weights
  opt1 <- optim(
    par = rep(start_val, ncol(EM)),
    fn = objfn, gr = gradfn,
    X = EM,
    method = method,
    control = list(maxit = 300, trace = 2), ...
  )
  alpha <- opt1$par
  wt <- exp(EM %*% alpha)
  wt_rs <- (wt / sum(wt)) * nrow(EM)

  # append weights to data
  data$weights <- NA
  data$weights[!ind] <- wt

  data$scaled_weights <- NA
  data$scaled_weights[!ind] <- wt_rs

  if (is.numeric(centered_colnames)) centered_colnames <- names(data)[centered_colnames]

  # Output
  list(
    data = data,
    centered_colnames = centered_colnames,
    nr_missing = nr_missing,
    ess = sum(wt)^2 / sum(wt^2),
    opt = opt1
  )
}


#' Plot MAIC weights in a histogram with key statistics in legend
#'
#' Generates a plot given the individuals weights with key summary in top right legend that includes
#' median weight, effective sample size (ESS), and reduction percentage (what percent ESS is reduced from the
#' original sample size).
#' There are two options of weights provided in \code{\link{estimate_weights}}: unscaled or scaled.
#' Scaled weights are relative to the original unit weights of each individual.
#' In other words, a scaled weight greater than 1 means that an individual carries more weight in the
#' re-weighted population than the original data and a scaled weight less than 1 means that an individual carries
#' less weight in the re-weighted population than the original data.
#'
#' @param wt a numeric vector of individual MAIC weights (derived using \code{\link{estimate_weights}})
#' @param bin_col a string, color for the bins of histogram
#' @param vline_col a string, color for the vertical line in the histogram
#' @param main_title a character string, main title of the plot
#'
#' @return a plot of unscaled or scaled weights
#' @importFrom graphics hist
#' @export

plot_weights <- function(wt, bin_col = "#6ECEB2", vline_col = "#688CE8", main_title = "Unscaled Individual Weights") {
  # calculate sample size and exclude NA from wt
  nr_na <- sum(is.na(wt))
  n <- length(wt) - nr_na
  wt <- na.omit(wt)

  # calculate effective sample size (ESS) and reduction of original sample
  ess <- (sum(wt)^2) / sum(wt^2)
  ess_reduct <- round((1 - (ess / n)) * 100, 2)

  # prepare legend
  plot_legend <- c(
    paste0("Median = ", round(median(wt), 4)),
    paste0("ESS = ", round(ess, 2)),
    paste0("Reduction% = ", ess_reduct)
  )
  plot_lty <- c(2, NA, NA)

  if (nr_na > 0) {
    plot_legend <- c(plot_legend, paste0("#Missing Weights = ", nr_na))
    plot_lty <- c(plot_lty, NA)
  }

  # plot
  par(mgp = c(2.3, 0.5, 0), cex.axis = 0.9, cex.lab = 0.95, bty = "n")
  hist(wt, border = "white", col = bin_col, main = main_title, breaks = 20, yaxt = "n", xlab = "")
  axis(2, las = 1)
  abline(v = median(wt), lty = 2, col = vline_col, lwd = 2)
  legend("topright", bty = "n", lty = plot_lty, cex = 0.8, legend = plot_legend)
}

#' Check to see if weights are optimized correctly
#'
#' This function checks to see if the optimization is done properly by checking the covariate averages
#' before and after adjustment.
#'
#' @param optimized object returned after calculating weights using \code{\link{estimate_weights}}
#' @param processed_agd a data frame, object returned after using \code{\link{process_agd}} or
#' aggregated data following the same naming convention
#' @param mean_digits number of digits for rounding mean columns in the output
#' @param prop_digits number of digits for rounding proportion columns in the output
#' @param sd_digits number of digits for rounding mean columns in the output
#'
#' @import DescTools
#'
#' @return data.frame of weighted and unweighted covariate averages of the IPD, and average of aggregate data
#' @export
#'
#' @examples
#' adsl <- read.csv(system.file("extdata", "adsl.csv",
#'   package = "maicplus",
#'   mustWork = TRUE
#' ))
#' adrs <- read.csv(system.file("extdata", "adrs.csv",
#'   package = "maicplus",
#'   mustWork = TRUE
#' ))
#' adtte <- read.csv(system.file("extdata", "adtte.csv",
#'   package = "maicplus",
#'   mustWork = TRUE
#' ))
#'
#' ### AgD
#' # Baseline aggregate data for the comparator population
#' target_pop <- read.csv(system.file("extdata", "aggregate_data_example_1.csv",
#'   package = "maicplus", mustWork = TRUE
#' ))
#' # target_pop2 <- read.csv(system.file("extdata", "aggregate_data_example_2.csv",
#' #                                     package = "maicplus", mustWork = TRUE))
#' # target_pop3 <- read.csv(system.file("extdata", "aggregate_data_example_3.csv",
#' #                                     package = "maicplus", mustWork = TRUE))
#'
#' # for time-to-event endpoints, pseudo IPD from digitalized KM
#' pseudo_ipd <- read.csv(system.file("extdata", "psuedo_IPD.csv",
#'   package = "maicplus",
#'   mustWork = TRUE
#' ))
#'
#' #### prepare data ----------------------------------------------------------
#' target_pop <- process_agd(target_pop)
#' # target_pop2 <- process_agd(target_pop2) # demo of process_agd in different scenarios
#' # target_pop3 <- process_agd(target_pop3) # demo of process_agd in different scenarios
#' adsl <- dummize_ipd(adsl, dummize_cols = c("SEX"), dummize_ref_level = c("Female"))
#' use_adsl <- center_ipd(ipd = adsl, agd = target_pop)
#'
#' match_res <- estimate_weights(
#'   data = use_adsl,
#'   centered_colnames = grep("_CENTERED$", names(use_adsl)),
#'   start_val = 0,
#'   method = "BFGS"
#' )
#'
#' check <- check_weights(
#'   optimized = match_res,
#'   processed_agd = target_pop
#' )
#'
#' print(check)

check_weights <- function(optimized, processed_agd) {
  ipd_with_weights <- optimized$data
  match_cov <- optimized$centered_colnames

  # if algorithm is correct, all centered columns should have a weighted summation to a very small number around zero
  num_check <- ipd_with_weights$weights %*% as.matrix(ipd_with_weights[, match_cov, drop = FALSE])
  num_check <- round(num_check, 4)

  # for reporting
  outdata <- data.frame(
    covariate = gsub("_CENTERED$", "", match_cov),
    match_stat = NA,
    internal_trial = NA,
    internal_trial_after_weighted = NA,
    external_trial = NA,
    sum_centered_IPD_with_weights = as.vector(num_check)
  )
  attr(outdata, "footer") <- list()
  ind_mean <- lapply(outdata$covariate, grep, x = names(processed_agd), value = TRUE) # find item that was matched by mean
  ind_mean <- sapply(ind_mean, function(ii) any(grepl("_MEAN$", ii)))
  outdata$match_stat <- ifelse(grepl("_MEDIAN$", outdata$covariate), "Median",
    ifelse(grepl("_SQUARED$", outdata$covariate), "SD",
      ifelse(ind_mean, "Mean", "Prop")
    )
  )
  outdata$covariate <- gsub("_MEDIAN|_SQUARED", "", outdata$covariate)
  # fill in corresponding agd data
  outdata$external_trial <- unlist(processed_agd[paste(outdata$covariate, toupper(outdata$match_stat), sep = "_")])

  # fill in stat from unweighted and weighted IPD
  for (ii in 1:nrow(outdata)) {
    covname <- outdata$covariate[ii]
    if (outdata$match_stat[ii] %in% c("Mean", "Prop")) {
      outdata$internal_trial[ii] <- mean(ipd_with_weights[[covname]], na.rm = TRUE)
      outdata$internal_trial_after_weighted[ii] <- weighted.mean(ipd_with_weights[[covname]], w = ipd_with_weights$weights, na.rm = TRUE)
    } else if (outdata$match_stat[ii] == "Median") {
      outdata$internal_trial[ii] <- quantile(ipd_with_weights[[covname]],
        probs = 0.5,
        na.rm = TRUE,
        type = 2,
        names = FALSE
      ) # SAS default
      outdata$internal_trial_after_weighted[ii] <- DescTools::Quantile(ipd_with_weights[[covname]],
        weights = ipd_with_weights$weights,
        probs = 0.5,
        na.rm = TRUE,
        type = 5,
        names = FALSE
      )
      # no IPD equals to reported AgD median
      msg_ind <- !any(ipd_with_weights[[covname]] == outdata$external_trial[ii], na.rm = TRUE)
      if (msg_ind) {
        msg_txt <- paste0(
          "For covariate ", covname, ", it was matched to AgD median, but there is no IPD identical to AgD median,",
          "hence median after weighted will not equal to AgD median exactly."
        )
        attr(outdata, "footer") <- c(attr(outdata, "footer"), msg_txt)
      }
    } else if (outdata$match_stat[ii] == "SD") {
      outdata$internal_trial[ii] <- sd(ipd_with_weights[[covname]], na.rm = TRUE)
      wm_squared <- weighted.mean(ipd_with_weights[[covname]]^2, w = ipd_with_weights$weights, na.rm = TRUE)
      ms_agd <- processed_agd[[paste0(outdata$covariate[ii], "_MEAN")]]^2
      outdata$internal_trial_after_weighted[ii] <- sqrt(wm_squared - ms_agd)
    }
  }

  # output
  class(outdata) <- c("maicplus_check_weights", "data.frame")
  outdata
}



#' Print method for Check Weights objects
#'
#' @param x object from [check_weights]
#' @param mean_digits number of digits for rounding mean columns in the output
#' @param prop_digits number of digits for rounding proportion columns in the output
#' @param sd_digits number of digits for rounding mean columns in the output
#' @param digits minimal number of significant digits, see [print.default].
#' @param ... further arguments to [print.data.frame]
#'
#' @describeIn check_weights Print method for check_weights objects
#' @export

print.maicplus_check_weights <- function(x, mean_digits = 2, prop_digits = 2, sd_digits = 3, digits = getOption("digits"), ...) {
  round_digits <- c("Mean" = mean_digits, "Prop" = prop_digits, "SD" = sd_digits)[x$match_stat]
  round_digits[is.na(round_digits)] <- digits

  x$external_trial <- round(x$external_trial, round_digits)
  x$internal_trial <- round(x$internal_trial, round_digits)
  x$internal_trial_after_weighted <- round(x$internal_trial_after_weighted, round_digits)

  print.data.frame(x, ...)
  footer <- unlist(attr(x, "footer"))
  if (length(footer)) {
    cat("\n")
    for (f in seq_along(footer)) {
      cat(paste0("[", f, "] ", footer[f]))
    }
  }
}
