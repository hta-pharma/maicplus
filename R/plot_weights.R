#' Plot MAIC weights in a histogram with key statistics in legend
#'
#' @param wt a numeric vector of individual MAIC weights (derived use in \code{\link{cal_weights}})
#' @param main.title a character string, main title of the plot
#'
#' @return a plot
#' @export
#'
#' @examples
plot_weights <- function(wt, main.title = "Unscaled Individual Weigths") {

  # calculate sample size and exclude NA from wt
  n <- length(wt)
  n.na <- sum(is.na(wt))
  wt <- na.omit(wt)

  # calculate effective sample size (ESS) and reduction of original sample
  ess <- (sum(wt)^2) / sum(wt^2)
  ess_reduct <- round((1 - (ess / n)) * 100, 2)

  # prepare legend
  plot.legend <- c(
    paste0("Median = ", round(median(wt), 4)),
    paste0("ESS = ", round(ess, 2)),
    paste0("Reduction% = ", ess_reduct)
  )
  plot.lty <- c(2, NA, NA)
  
  if (n.na > 0){
    plot.legend <- c(plot.legend, paste0("#Missing Weights = ", n.na))
    plot.lty <- c(plot.lty, NA)
  }

  # plot
  par(mfrow = c(1, 1), family = "HersheySerif", mgp = c(2.3, 0.5, 0), cex.axis = 0.9, cex.lab = 0.95, bty = "n")
  hist(wt, border = "white", col = "#6ECEB2", main = main.title, breaks = 20, yaxt = "n", xlab = "")
  axis(2, las = 1)
  abline(v = median(wt), lty = 2, col = "#688CE8", lwd = 2)
  legend("topright", bty = "n", lty = plot.lty, cex = 0.8, legend = plot.legend)
}
