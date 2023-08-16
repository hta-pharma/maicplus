# Functions for matching step: estimation of individual weights

#' Derive individual weights in the matching step of MAIC
#'
#' This function takes individual patient data (IPD) with centered covariates
#' (effect modifiers and/or prognostic variables) as input and generates
#' weights for each individual in IPD trial to match the covariates in aggregate data.
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
  outdata <- list(
    data = data,
    centered_colnames = centered_colnames,
    nr_missing = nr_missing,
    ess = sum(wt)^2 / sum(wt^2),
    opt = opt1
  )

  class(outdata) <- c("maicplus_estimate_weights", "list")
  outdata
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
#' @export

plot_weights_base <- function(weighted_data,
                              bin_col, vline_col, main_title,
                              scaled_weights) {
  if (scaled_weights) {
    wt <- weighted_data$data$scaled_weights
  } else {
    wt <- weighted_data$data$weights
  }

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

#' Plot MAIC weights in a histogram with key statistics in legend using ggplot
#'
#' Generates a ggplot histogram of weights. Default is to plot both unscaled and scaled weights on a same graph.
#'
#' @param weighted_data object returned after calculating weights using [estimate_weights]
#' @param bin_col a string, color for the bins of histogram
#' @param vline_col a string, color for the vertical line in the histogram
#' @param main_title Name of scaled weights plot and unscaled weights plot, respectively.
#' @param bins number of bin parameter to use
#' @param print_caption print a footnote message related to ESS from the NICE survey 2021
#' @param caption_width width that is passed onto str_wrap function
#'
#' @return a plot of unscaled and scaled weights
#' @export

plot_weights_ggplot <- function(weighted_data, bin_col, vline_col,
                                main_title, bins, print_caption, caption_width) {
  # check if survminer package is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is needed to run this function")
  }

  wt_data0 <- weighted_data$data[, c("weights", "scaled_weights")]
  colnames(wt_data0) <- c(main_title[1], main_title[2])
  wt_data <- stack(wt_data0)
  wt_data$median <- ifelse(wt_data$ind == main_title[1], median(wt_data0[, main_title[1]]), median(wt_data0[, main_title[2]]))

  summ <- aggregate(wt_data$values, list(wt_data$ind), median)
  colnames(summ) <- c("ind", "lab")
  summ$lab <- paste0(
    "Median = ", round(summ$lab, 4),
    "\nESS = ", round(weighted_data$ess, 2),
    "\nReduction% = ", round((1 - (weighted_data$ess / dim(weighted_data$data)[1])) * 100, 2)
  )

  hist_plot <- ggplot2::ggplot(wt_data) +
    ggplot2::geom_histogram(ggplot2::aes(values), bins = bins, color = bin_col, fill = bin_col) +
    ggplot2::geom_vline(aes(xintercept = median),
      color = vline_col,
      linetype = "dashed", linewidth = 0.5
    ) +
    theme_bw() +
    ggplot2::facet_wrap(~ind, ncol = 1) + # gives the two plots (one on top of the other)
    ggplot2::geom_text(data = summ, aes(label = lab), x = Inf, y = Inf, hjust = 1, vjust = 1, size = 3) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 16)
    ) +
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("Weight")

  if (print_caption == TRUE) {
    print_text <- "In most applications, weighting considerably reduces the effective sample size from the original AC sample size. The median percentage reduction is 58% (range: 7.9%–94.1%; interquartile range: 42.2%–74.2%). The final effective sample sizes are also representative of those in the technology appraisals, which are also small (median: 80; range: 4.8–639; interquartile range: 37–174). Therefore, an ESS reduction up to ~60% is not unexpected based on the 2021 survey, whereas a reduction of >75% is less common and it may be considered suboptimal."
    hist_plot <- hist_plot + ggplot2::labs(caption = stringr::str_wrap(print_text, width = caption_width)) +
      ggplot2::theme(
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0)
      )
  }
  return(hist_plot)
}


#' Plot method for Estimate Weights objects
#'
#' The plot function displays individuals weights with key summary in top right legend that includes
#' median weight, effective sample size (ESS), and reduction percentage (what percent ESS is reduced from the
#' original sample size). There are two options of plotting: base R plot and ggplot. The default
#' for base R plot is to plot unscaled and scaled separately. The default
#' for ggplot is to plot unscaled and scaled weights on a same plot.
#'
#' @param x object from [estimate_weights]
#' @param ggplot indicator to print base weights plot or ggplot weights plot
#' @param bin_col a string, color for the bins of histogram
#' @param vline_col a string, color for the vertical line in the histogram
#' @param main_title title of the plot. For ggplot, name of scaled weights plot and unscaled weights plot, respectively.
#' @param scaled_weights (base plot only) an indicator for using scaled weights instead of regular weights
#' @param bins (ggplot only) number of bin parameter to use
#' @param print_caption (ggplot only) print a footnote message related to ESS from the NICE survey 2021
#' @param caption_width (ggplot only) width that is passed onto str_wrap function
#'
#' @examples
#' @describeIn estimate_weights Plot method for estimate_weights objects
#' @export

plot.maicplus_estimate_weights <- function(x, ggplot = FALSE,
                                           bin_col = NULL, vline_col = NULL,
                                           main_title = NULL, scaled_weights = TRUE,
                                           bins = 50, print_caption = FALSE,
                                           caption_width = 80) {
  if (!ggplot) {
    if (is.null(main_title)) main_title <- "Scaled Individual Weights"
    if (is.null(bin_col)) bin_col <- "#6ECEB2"
    if (is.null(vline_col)) vline_col <- "#688CE8"
  } else {
    if (is.null(main_title)) main_title <- c("Scaled Individual Weights", "Unscaled Individual Weights")
    if (is.null(bin_col)) bin_col <- "black"
    if (is.null(vline_col)) vline_col <- "red"
  }

  if (!ggplot) {
    plot_weights_base(x, bin_col, vline_col, main_title, scaled_weights)
  } else {
    plot_weights_ggplot(x, bin_col, vline_col, main_title, bins, print_caption, caption_width)
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
#'   weighted_data = match_res,
#'   processed_agd = target_pop
#' )
#'
#' print(check)
#'
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

#' Prints Note on Expected Sample Size Reduction
#'
#' @param width Number of characters to break string into new lines (`\n`).
#'
#' @return A character string
#' @export
#'
#' @examples
#' ess_footnote_text(width = 80)
ess_footnote_text <- function(width = 0.9 * getOption("width")) {
  text <- "An ESS reduction up to ~60% is not unexpected based on the 2021 survey of NICE's technology appraisals
(https://onlinelibrary.wiley.com/doi/full/10.1002/jrsm.1511), whereas a reduction of >75% is less common
and it may be considered sub optimal."
  paste0(strwrap(text, width = width), collapse = "\n")
}