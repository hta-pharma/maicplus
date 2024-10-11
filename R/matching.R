# Functions for matching step: estimation of individual weights

# functions to be exported ---------------------------------------

#' Derive individual weights in the matching step of MAIC
#'
#' Assuming data is properly processed, this function takes individual patient data (IPD) with centered covariates
#' (effect modifiers and/or prognostic variables) as input, and generates weights for each individual in IPD trial to
#' match the covariates in aggregate data.
#'
#' @param data a numeric matrix, centered covariates of IPD, no missing value in any cell is allowed
#' @param centered_colnames a character or numeric vector (column indicators) of centered covariates
#' @param start_val a scalar, the starting value for all coefficients of the propensity score regression
#' @param method a string, name of the optimization algorithm (see 'method' argument of \code{base::optim()}) The
#'   default is `"BFGS"`, other options are `"Nelder-Mead"`, `"CG"`, `"L-BFGS-B"`, `"SANN"`, and `"Brent"`
#' @param n_boot_iteration an integer, number of bootstrap iterations. By default is NULL which means bootstrapping
#'   procedure will not be triggered, and hence the element `"boot"` of output list object will be NULL.
#' @param set_seed_boot a scalar, the random seed for conducting the bootstrapping, only relevant if
#'   \code{n_boot_iteration} is not NULL. By default, use seed 1234
#' @param boot_strata a character vector of column names in \code{data} that defines the strata for bootstrapping.
#'   This ensures that samples are drawn proportionally from each defined stratum. If \code{NULL},
#'   no stratification during bootstrapping process. By default, it is "ARM"
#' @param ... Additional `control` parameters passed to [stats::optim].
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
#'   \item{boot_strata}{'strata' from a boot::boot object}
#'   \item{boot_seed}{column names in \code{data} of the stratification factors}
#'   \item{boot}{a n by 2 by k array or NA, where n equals to number of rows in \code{data}, and k equals
#'      \code{n_boot_iteration}. The 2 columns in the second dimension include a column of numeric indexes of the rows
#'      in \code{data} that are selected at a bootstrapping iteration and a column of weights. \code{boot} is NA when
#'      argument \code{n_boot_iteration} is set as NULL
#'   }
#' }
#' @importFrom boot boot
#' @examples
#' data(centered_ipd_sat)
#' centered_colnames <- grep("_CENTERED", colnames(centered_ipd_sat), value = TRUE)
#' weighted_data <- estimate_weights(data = centered_ipd_sat, centered_colnames = centered_colnames)
#' \donttest{
#' # To later estimate bootstrap confidence intervals, we calculate the weights
#' # for the bootstrap samples:
#' weighted_data_boot <- estimate_weights(
#'   data = centered_ipd_sat, centered_colnames = centered_colnames, n_boot_iteration = 100
#' )
#' }
#' @export

estimate_weights <- function(data,
                             centered_colnames = NULL,
                             start_val = 0,
                             method = "BFGS",
                             n_boot_iteration = NULL,
                             set_seed_boot = 1234,
                             boot_strata = "ARM",
                             ...) {
  # pre check
  ch1 <- is.data.frame(data)
  if (!ch1) {
    stop("'data' is not a data.frame")
  }

  ch2 <- (!is.null(centered_colnames))
  if (ch2 && is.numeric(centered_colnames)) {
    ch2b <- any(centered_colnames < 1 | centered_colnames > ncol(data))
    if (ch2b) {
      stop("specified centered_colnames are out of bound")
    }
  } else if (ch2 && is.character(centered_colnames)) {
    ch2b <- !all(centered_colnames %in% names(data))
    if (ch2b) {
      stop("1 or more specified centered_colnames are not found in 'data'")
    }
  } else {
    stop("'centered_colnames' should be either a numeric or character vector")
  }

  ch3 <- sapply(centered_colnames, function(ii) {
    !is.numeric(data[[ii]])
  })
  if (any(ch3)) {
    stop(paste0(
      "following columns of 'data' are not numeric for the calculation:",
      paste(which(ch3), collapse = ",")
    ))
  }

  if (!is.null(boot_strata)) {
    ch4 <- boot_strata %in% names(data)
    if (!all(ch4)) {
      stop("Some variables in boot_strata are not in data: ", toString(boot_strata[!ch4]))
    }
  }

  # prepare data for optimization
  if (is.null(centered_colnames)) centered_colnames <- seq_len(ncol(data))
  EM <- data[, centered_colnames, drop = FALSE]
  ind <- apply(EM, 1, function(xx) any(is.na(xx)))
  nr_missing <- sum(ind)
  rows_with_missing <- which(ind)
  EM <- as.matrix(EM[!ind, , drop = FALSE])

  # estimate weights
  opt1 <- optimise_weights(matrix = EM, par = rep(start_val, ncol(EM)), method = method, ...)
  alpha <- opt1$alpha
  wt <- opt1$wt
  wt_rs <- (wt / sum(wt)) * nrow(EM)

  # bootstrapping
  outboot <- if (is.null(n_boot_iteration)) {
    boot_seed <- NULL
    boot_strata_out <- NULL
    NULL
  } else {
    # Make sure to leave '.Random.seed' as-is on exit
    old_seed <- globalenv()$.Random.seed
    on.exit(suspendInterrupts(set_random_seed(old_seed)))
    set.seed(set_seed_boot)

    if (!is.null(boot_strata)) {
      use_strata <- interaction(data[!ind, boot_strata])
    } else {
      use_strata <- rep(1, nrow(EM))
    }
    boot_statistic <- function(d, w) optimise_weights(d[w, ], par = alpha, method = method, ...)$wt[, 1]
    boot_out <- boot::boot(EM, statistic = boot_statistic, R = n_boot_iteration, strata = use_strata)

    boot_array <- array(dim = list(nrow(EM), 2, n_boot_iteration))
    dimnames(boot_array) <- list(sampled_patient = NULL, c("rowid", "weight"), bootstrap_iteration = NULL)
    boot_array[, 1, ] <- t(boot.array(boot_out, TRUE))
    boot_array[, 2, ] <- t(boot_out$t)
    boot_seed <- boot_out$seed
    boot_strata_out <- boot_out$strata
    boot_array
  }

  # append weights to data
  data$weights <- NA
  data$weights[!ind] <- wt

  data$scaled_weights <- NA
  data$scaled_weights[!ind] <- wt_rs

  if (is.numeric(centered_colnames)) centered_colnames <- names(data)[centered_colnames]

  # Output
  outdata <- list(
    data = data,
    centered_colnames = centered_colnames,
    nr_missing = nr_missing,
    ess = sum(wt)^2 / sum(wt^2),
    opt = opt1$opt,
    boot = outboot,
    boot_seed = boot_seed,
    boot_strata = boot_strata_out,
    rows_with_missing = rows_with_missing
  )

  class(outdata) <- c("maicplus_estimate_weights", "list")
  outdata
}


#' Estimate weights using `optim`
#'
#' @param matrix Matrix of data to be used for estimating weights
#' @param par Vector of starting values for the parameters with length equal to the number of columns in `matrix`
#' @param method Method parameter passed to [stats::optim]
#' @param ... Additional `control` parameters passed to [stats::optim]
#'
#' @return List containing estimated `alpha` values and `wt` weights for all rows of matrix
#' @noRd
optimise_weights <- function(matrix,
                             par = rep(0, ncol(matrix)),
                             method = "BFGS",
                             maxit = 300,
                             trace = 0,
                             ...) {
  if (!all(is.numeric(par) || is.finite(par), length(par) == ncol(matrix))) {
    stop("par must be a numeric vector with finite values of length equal to the number of columns in 'matrix'")
  }
  opt1 <- optim(
    par = par,
    fn = function(alpha, X) sum(exp(X %*% alpha)),
    gr = function(alpha, X) colSums(sweep(X, 1, exp(X %*% alpha), "*")),
    X = matrix,
    method = method,
    control = list(maxit = maxit, trace = trace, ...)
  )
  if (opt1$convergence != 0) {
    warning(
      "optim() did not converge. ",
      opt1$message,
      "\nSee ?optim for more information on convergence code: ", opt1$convergence
    )
  }
  list(
    opt = opt1,
    alpha = opt1$par,
    wt = exp(matrix %*% opt1$par)
  )
}

#' Calculate Statistics for Weight Plot Legend
#'
#' Calculates ESS reduction and median weights which is used to create legend for weights plot
#'
#' @param weighted_data object returned after calculating weights using [estimate_weights]
#'
#' @return list of ESS, ESS reduction, median value of scaled and unscaled weights, and missing count
#' @examples
#' data("weighted_sat")
#' calculate_weights_legend(weighted_sat)
#' @export
#' @keywords internal

calculate_weights_legend <- function(weighted_data) {
  if (!inherits(weighted_data, "maicplus_estimate_weights")) {
    stop("weighted_data must be class `maicplus_estimate_weights` generated by estimate_weights()")
  }
  ess <- weighted_data$ess
  wt <- weighted_data$data$weights
  wt_scaled <- weighted_data$data$scaled_weights

  # calculate sample size and exclude NA from wt
  nr_na <- sum(is.na(wt))
  n <- length(wt) - nr_na
  wt <- na.omit(wt)
  wt_scaled <- na.omit(wt_scaled)

  # calculate ess reduction and median weights
  ess_reduction <- (1 - (ess / n)) * 100
  wt_median <- median(wt)
  wt_scaled_median <- median(wt_scaled)

  list(
    ess = round(ess, 2),
    ess_reduction = round(ess_reduction, 2),
    wt_median = round(wt_median, 4),
    wt_scaled_median = round(wt_scaled_median, 4),
    nr_na = nr_na
  )
}

#' Plot MAIC weights in a histogram with key statistics in legend
#'
#' Generates a base R histogram of weights. Default is to plot either unscaled or scaled weights and not both.
#'
#' @param weighted_data object returned after calculating weights using [estimate_weights]
#' @param bin_col a string, color for the bins of histogram
#' @param vline_col a string, color for the vertical line in the histogram
#' @param main_title title of the plot
#' @param scaled_weights an indicator for using scaled weights instead of regular weights
#'
#' @return a plot of unscaled or scaled weights
#' @examples
#' plot_weights_base(weighted_sat,
#'   bin_col = "#6ECEB2",
#'   vline_col = "#688CE8",
#'   main_title = c("Scaled Individual Weights", "Unscaled Individual Weights"),
#'   scaled_weights = TRUE
#' )
#' @export

plot_weights_base <- function(weighted_data,
                              bin_col, vline_col, main_title,
                              scaled_weights) {
  weights_stat <- calculate_weights_legend(weighted_data)

  if (scaled_weights) {
    wt <- weighted_data$data$scaled_weights
    median_wt <- weights_stat$wt_scaled_median
  } else {
    wt <- weighted_data$data$weights
    median_wt <- weights_stat$wt_median
  }

  # prepare legend
  plot_legend <- c(
    paste0("Median = ", median_wt),
    paste0("ESS = ", weights_stat$ess),
    paste0("Reduction% = ", weights_stat$ess_reduct)
  )
  plot_lty <- c(2, NA, NA)

  if (weights_stat$nr_na > 0) {
    plot_legend <- c(plot_legend, paste0("#Missing Weights = ", weights_stat$nr_na))
    plot_lty <- c(plot_lty, NA)
  }

  # plot
  original_par <- par(mgp = c(2.3, 0.5, 0), cex.axis = 0.9, cex.lab = 0.95, bty = "n")
  on.exit(par(original_par))
  hist(wt, border = "white", col = bin_col, main = main_title, breaks = 20, yaxt = "n", xlab = "")
  axis(2, las = 1)
  abline(v = median(wt), lty = 2, col = vline_col, lwd = 2)
  legend("topright", bty = "n", lty = plot_lty, cex = 0.8, legend = plot_legend)
}

#' Plot MAIC weights in a histogram with key statistics in legend using `ggplot2`
#'
#' Generates a `ggplot` histogram of weights. Default is to plot both unscaled and scaled weights on a same graph.
#'
#' @param weighted_data object returned after calculating weights using [estimate_weights]
#' @param bin_col a string, color for the bins of histogram
#' @param vline_col a string, color for the vertical line in the histogram
#' @param main_title Name of scaled weights plot and unscaled weights plot, respectively.
#' @param bins number of bin parameter to use
#'
#' @return a plot of unscaled and scaled weights
#' @examples
#' if (requireNamespace("ggplot2")) {
#'   plot_weights_ggplot(weighted_sat,
#'     bin_col = "#6ECEB2",
#'     vline_col = "#688CE8",
#'     main_title = c("Scaled Individual Weights", "Unscaled Individual Weights"),
#'     bins = 50
#'   )
#' }
#' @export

plot_weights_ggplot <- function(weighted_data, bin_col, vline_col,
                                main_title,
                                bins) {
  # check if ggplot2 package is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is needed to run this function")
  }

  weights_stat <- calculate_weights_legend(weighted_data)

  # prepare dataset to use in ggplot
  wt_data0 <- weighted_data$data[, c("weights", "scaled_weights")]
  colnames(wt_data0) <- main_title
  wt_data <- stack(wt_data0)
  wt_data$median <- ifelse(wt_data$ind == main_title[1],
    weights_stat$wt_median, weights_stat$wt_scaled_median
  )

  # create legend data
  lab <- with(weights_stat, {
    lab <- c(paste0("Median = ", wt_median), paste0("Median = ", wt_scaled_median))
    lab <- paste0(lab, "\nESS = ", ess, "\nReduction% = ", ess_reduction)
    if (nr_na > 0) lab <- paste0(lab, "\n#Missing Weights = ", nr_na)
    lab
  })
  legend_data <- data.frame(ind = main_title, lab = lab)

  hist_plot <- ggplot2::ggplot(wt_data) +
    ggplot2::geom_histogram(ggplot2::aes_string(x = "values"), bins = bins, color = bin_col, fill = bin_col) +
    ggplot2::geom_vline(ggplot2::aes_string(xintercept = "median"),
      color = vline_col,
      linetype = "dashed"
    ) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~ind, ncol = 1) +
    ggplot2::geom_text(
      data = legend_data,
      ggplot2::aes_string(label = "lab"), x = Inf, y = Inf, hjust = 1, vjust = 1, size = 3
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("Weight")

  return(hist_plot)
}


#' Plot method for Estimate Weights objects
#'
#' The plot function displays individuals weights with key summary in top right legend that includes
#' median weight, effective sample size (ESS), and reduction percentage (what percent ESS is reduced from the
#' original sample size). There are two options of plotting: base R plot and `ggplot`. The default
#' for base R plot is to plot unscaled and scaled separately. The default
#' for `ggplot` is to plot unscaled and scaled weights on a same plot.
#'
#' @param x object from [estimate_weights]
#' @param ggplot indicator to print base weights plot or `ggplot` weights plot
#' @param bin_col a string, color for the bins of histogram
#' @param vline_col a string, color for the vertical line in the histogram
#' @param main_title title of the plot. For ggplot, name of scaled weights plot and unscaled weights plot, respectively.
#' @param scaled_weights (base plot only) an indicator for using scaled weights instead of regular weights
#' @param bins (`ggplot` only) number of bin parameter to use
#'
#' @examples
#' plot(weighted_sat)
#'
#' if (requireNamespace("ggplot2")) {
#'   plot(weighted_sat, ggplot = TRUE)
#' }
#' @describeIn estimate_weights Plot method for estimate_weights objects
#' @export

plot.maicplus_estimate_weights <- function(x, ggplot = FALSE,
                                           bin_col = "#6ECEB2", vline_col = "#688CE8",
                                           main_title = NULL,
                                           scaled_weights = TRUE,
                                           bins = 50, ...) {
  if (ggplot) {
    if (is.null(main_title)) main_title <- c("Scaled Individual Weights", "Unscaled Individual Weights")
    plot_weights_ggplot(x, bin_col, vline_col, main_title, bins)
  } else {
    if (is.null(main_title)) {
      main_title <- ifelse(scaled_weights, "Scaled Individual Weights", "Unscaled Individual Weights")
    }
    plot_weights_base(x, bin_col, vline_col, main_title, scaled_weights)
  }
}


#' Check to see if weights are optimized correctly
#'
#' This function checks to see if the optimization is done properly by checking the covariate averages
#' before and after adjustment.
#'
#' @param weighted_data object returned after calculating weights using \code{\link{estimate_weights}}
#' @param processed_agd a data frame, object returned after using \code{\link{process_agd}} or
#' aggregated data following the same naming convention
#'
#' @examples
#' data(weighted_sat)
#' data(agd)
#' check_weights(weighted_sat, process_agd(agd))
#'
#' @import DescTools
#'
#' @return data.frame of weighted and unweighted covariate averages of the IPD,
#' average of aggregate data, and sum of inner products of covariate \eqn{x_i} and the weights (\eqn{exp(x_i\beta)})
#' @export

check_weights <- function(weighted_data, processed_agd) {
  ipd_with_weights <- weighted_data$data
  match_cov <- weighted_data$centered_colnames

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
  # find item that was matched by mean
  ind_mean <- lapply(outdata$covariate, grep, x = names(processed_agd), value = TRUE)
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
  for (ii in seq_len(nrow(outdata))) {
    covname <- outdata$covariate[ii]
    if (outdata$match_stat[ii] %in% c("Mean", "Prop")) {
      outdata$internal_trial[ii] <- mean(ipd_with_weights[[covname]], na.rm = TRUE)
      outdata$internal_trial_after_weighted[ii] <- weighted.mean(
        ipd_with_weights[[covname]],
        w = ipd_with_weights$weights, na.rm = TRUE
      )
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
#' @describeIn check_weights Print method for check_weights objects
#' @export

print.maicplus_check_weights <- function(x,
                                         mean_digits = 2,
                                         prop_digits = 2,
                                         sd_digits = 3,
                                         digits = getOption("digits"), ...) {
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

#' Note on Expected Sample Size Reduction
#'
#' @param width Number of characters to break string into new lines (`\n`).
#'
#' @return A character string
#' @keywords internal
ess_footnote_text <- function(width = 0.9 * getOption("width")) {
  text <- "An ESS reduction up to ~60% is not unexpected based on the 2021 survey of NICE's technology appraisals
(https://onlinelibrary.wiley.com/doi/full/10.1002/jrsm.1511), whereas a reduction of >75% is less common
and it may be considered suboptimal."
  paste0(strwrap(text, width = width), collapse = "\n")
}
