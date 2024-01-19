#' helper function: sort out a nice report table to summarize survival analysis results
#' for `maic_tte_unanchor`
#'
#' @param coxobj returned object from \code{\link[survival]{coxph}}
#' @param medSurvobj returned object from \code{\link{medSurv_makeup}}
#' @param tag a string, by default NULL, if specified, an extra 1st column is created in the output
#'
#' @return a data frame with sample size, incidence rate, median survival time with 95% CI, hazard ratio estimate with
#' 95% CI and Wald test of hazard ratio
#'
#' @export

report_table <- function(coxobj, medSurvobj, tag = NULL) {
  hr_res <- format(round(summary(coxobj)$conf.int[-2], 2), nsmall = 2)
  hr_res <- paste0(hr_res[1], "[", hr_res[2], ";", hr_res[3], "]")
  hr_pval <- format(round(summary(coxobj)$waldtest[3], 3), nsmall = 3)
  if (hr_pval == "0.000") hr_pval <- "<0.001"

  meds_report <- format(round(medSurvobj[, c("median", "0.95LCL", "0.95UCL")], 1), nsmall = 1)
  meds_report <- apply(meds_report, 1, function(xx) paste0(xx[1], "[", xx[2], ";", xx[3], "]"))

  desc_res <- cbind(medSurvobj[, "treatment", drop = FALSE],
    data.frame(N = round(medSurvobj$n.max, 1)),
    "n.events(%)" = paste0(
      round(medSurvobj$events, 1), "(",
      format(round(medSurvobj$events * 100 / medSurvobj$n.max, 1), nsmall = 1), ")"
    ),
    "median[95% CI]" = meds_report
  )

  desc_res <- cbind(desc_res[c(2, 1), ], "HR[95% CI]" = c(hr_res, ""), "p-Value" = c(hr_pval, ""))

  if (!is.null(tag)) {
    desc_res <- cbind(data.frame(Matching = rep(tag, nrow(desc_res))), desc_res)
    desc_res$Matching[duplicated(desc_res$Matching)] <- ""
  }
  desc_res
}
