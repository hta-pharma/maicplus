#' Basic Kaplan Meier (KM) plot function
#'
#' This function can generate a basic KM plot with or without risk set table appended at the bottom.
#' In a single plot, it can include up to 3 KM curves. This depends on number of levels in 'treatment' column in the input data.frame \code{kmdat}
#'
#' @param kmdat a data.frame, must consist 'treatment', 'time' (unit in days), 'n.risk', 'censor', 'surv', similar to an output from \code{maicplus:::survfit_makeup}
#' @param time_scale a string, time unit of median survival time, taking a value of 'year', 'month', 'week' or 'day'
#' @param show_risk_set logical, show risk set table or not, TRUE by default
#' @param main_title a string, main title of the KM plot
#' @param subplot_heights a numeric vector, heights argument to \code{graphic::layout()},NULL by default which means user will use the default setting
#' @param suppress_plot_layout logical, suppress the layout setting in this function so that user can specify layout outside of the function, FALSE by default
#' @param use_colors a character vector of length up to 3, colors to the KM curves, it will be passed to 'col' of \code{lines()}
#' @param use_line_types a numeric vector of length up to 3, line type to the KM curves, it will be passed to 'lty' of \code{lines()}
#'
#' @example basic_kmplot_ex.R
#'
#' @return a KM plot with or without risk set table appended at the bottom, with up to 3 KM curves
#' @export

basic_kmplot <- function(kmdat, time_scale, show_risk_set = TRUE,
                         main_title = "Kaplan-Meier Curves",
                         subplot_heights = NULL,
                         suppress_plot_layout = FALSE,
                         use_colors = NULL,
                         use_line_types = NULL) {

  time_unit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)

  # precheck
  if (!length(subplot_heights) %in% c(0, (1 + show_risk_set))) stop(paste("length of subplot_heights should be", (1 + show_risk_set)))
  if (!is.factor(kmdat$treatment)) stop("kmdat$treatment needs to be a factor, its levels will be used in legend and title, first level is comparator")
  if (nlevels(kmdat$treatment) > 3) stop("kmdat$treatment cannot have more than 3 levels")
  if (is.null(time_grid) & show_risk_set) stop("please provide a numeric vector as time_grid to show risk set table")

  # set up x axis (time)
  if (is.null(time_grid)) {
    max_t <- max(kmdat$time)
    t_range <- c(0, (max_t / time_unit[[time_scale]]) * 1.07)
  } else {
    t_range <- c(0, max(time_grid))
  }

  # plat layout in par
  if (!suppress_plot_layout) {
    nr_subplot <- (1 + show_risk_set)
    if (is.null(subplot_heights)) subplot_heights <- c(7, 0.7 + nlevels(kmdat$treatment) * 0.7, 0.8)
    layout_mat <- matrix(1:nr_subplot, ncol = 1)
    layout(layout_mat, heights = subplot_heights)
  }

  # plot cosmetic setup
  if (is.null(use_line_types)) {
    use_lty <- c(1, 1, 2)
    use_lwd <- c(1.2, 2, 1.2)
  } else {
    use_lty <- use_line_types
    use_lwd <- c(1.2, 2, 1.2)
  }

  if (is.null(use_colors)) {
    use_col <- c("#5450E4", "#00857C", "#6ECEB2")
  } else {
    use_col <- use_colors
  }

  ## : first subplot: KM curve -------------------
  # base plot
  par(bty = "n", tcl = -0.15, mgp = c(1.8, 0.4, 0), cex.lab = 0.85, cex.axis = 0.8, cex.main = 0.9, mar = c(3, 4, 5, 1))
  plot(0, 0,
       type = "n", xlab = paste0("Time in ", time_scale), ylab = "Survival Probability",
       ylim = c(0, 1), xlim = t_range, yaxt = "n", xaxt = "n",
       main = main_title
  )
  axis(2, las = 1)
  if (!is.null(time_grid)) {
    axis(1, at = time_grid)
  } else {
    axis(1)
  }

  # add km line
  for (ii in 1:nlevels(kmdat$treatment)) {
    tmpkmdat <- kmdat[as.numeric(kmdat$treatment) == ii, , drop = FALSE]
    lines(
      y = tmpkmdat$surv,
      x = (tmpkmdat$time / time_unit[[time_scale]]),
      col = use_col[ii],
      lty = use_lty[ii],
      lwd = use_lwd[ii],
      type = "s"
    )
    tmpid <- tmpkmdat$censor == 1
    points(
      y = tmpkmdat$surv[tmpid],
      x = (tmpkmdat$time[tmpid] / time_unit[[time_scale]]),
      col = use_col[ii],
      pch = 3,
      cex = 0.65
    )
  }

  ## : second subplot: risk set table (if applicable) -------------------
  if (show_risk_set) {
    # add legend, with treatment index
    legend("topright",
           bty = "n",
           cex = 0.8,
           lty = use_lty[1:nlevels(kmdat$treatment)],
           lwd = use_lwd[1:nlevels(kmdat$treatment)],
           col = use_col[1:nlevels(kmdat$treatment)],
           legend = paste0("(T", 1:nlevels(kmdat$treatment), ") ", levels(kmdat$treatment))
    )

    # add risk set table
    par(bty = "n", tcl = -0.15, mgp = c(1.8, 0.4, 0), mar = c(1, 4, 0, 1))
    plot(0, 0,
         type = "n", xlab = "", ylab = "", main = NULL,
         ylim = c(nlevels(kmdat$treatment) + 1.2, -0.5),
         xlim = t_range,
         yaxt = "n", xaxt = "n"
    )
    axis(2,
         at = 1:nlevels(kmdat$treatment), labels = paste0("T", 1:nlevels(kmdat$treatment)),
         line = NA, lty = "blank", las = 1
    )

    for (ii in 1:nlevels(kmdat$treatment)) {
      tmpkmdat <- kmdat[as.numeric(kmdat$treatment) == ii, , drop = FALSE]
      tmptime <- (tmpkmdat$time / time_unit[[time_scale]])
      tmpnr <- sapply(time_grid, function(kk) {
        tmpid <- which(tmptime > kk)
        if (length(tmpid) == 0) {
          if (min(tmpkmdat$n.risk) == 0) {
            tmpid <- which.min(tmpkmdat$n.risk)[1]
          } else {
            tmpid <- NULL
          }
        }
        tout <- ifelse(is.null(tmpid), "n/a", tmpkmdat$n.risk[tmpid])
      })
      text(0, 0, labels = "Number at risk", pos = 4, cex = 0.8, offset = -0.8)
      text(
        y = rep(ii, length(time_grid)),
        x = time_grid,
        labels = tmpnr,
        col = use_col[ii],
        cex = 0.75
      )
      text(0, nlevels(kmdat$treatment) + 1,
           pos = 4, cex = 0.7, offset = -0.8, col = "gray30",
           labels = "Note: Number at risk for adjusted/weighted treament arm is the sum of individual weight at risk."
      )
    }
  } else {
    # add simple
    legend("topright",
           bty = "n",
           cex = 0.8,
           lty = use_lty[1:nlevels(kmdat$treatment)],
           lwd = use_lwd[1:nlevels(kmdat$treatment)],
           col = use_col[1:nlevels(kmdat$treatment)],
           legend = levels(kmdat$treatment)
    )
  }
}
