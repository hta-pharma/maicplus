#' helper function: sort out a nice report table to summarize survival analysis results
#'
#' @param coxobj returned object from \code{\link[survival]{coxph}}
#' @param medSurvobj returned object from \code{\link{medSurv_makeup}}
#' @param tag a string, by default NULL, if specified, an extra 1st column is created in the output
#'
#' @return a data frame with sample size, incidence rate, median survival time with 95% CI, hazard ratio estimate with
#' 95% CI and Wald test of hazard ratio
#'
#' @export

report_table_tte <- function(coxobj, medSurvobj, tag = NULL) {
  # descriptive part
  N <- medSurvobj$n.max
  N.EVNT <- medSurvobj$events
  N.EVNT.PERC <- round(N.EVNT * 100 / N, 1)
  N <- round(N, 1) # not necessary integer
  N.EVNT <- round(N.EVNT, 1) # not necessary integer
  meds_report <- format(round(medSurvobj[, c("median", "0.95LCL", "0.95UCL")], 1), nsmall = 1)
  meds_report <- apply(meds_report, 1, function(xx) paste0(xx[1], "[", xx[2], ";", xx[3], "]"))

  # inferential part
  hr_res <- format(round(summary(coxobj)$conf.int[-2], 2), nsmall = 2)
  hr_res <- paste0(hr_res[1], "[", hr_res[2], ";", hr_res[3], "]")
  hr_pval <- format(round(summary(coxobj)$waldtest[3], 3), nsmall = 3)
  if (hr_pval == "0.000") hr_pval <- "<0.001"

  # assemble the table
  desc_res <- cbind(
    medSurvobj[, "treatment", drop = FALSE],
    data.frame(
      "N" = N,
      "n.events(%)" = paste0(N.EVNT, "(", format(N.EVNT.PERC, nsmall = 1), ")"),
      "median[95% CI]" = meds_report
    )
  )
  desc_res <- cbind(desc_res[c(2, 1), ], "HR[95% CI]" = c(hr_res, ""), "p-Value" = c(hr_pval, ""))
  names(desc_res)[3] <- "n.events(%)" # cbind changes column name

  # add first tag column if applicable
  if (!is.null(tag)) {
    desc_res <- cbind(data.frame(Matching = rep(tag, nrow(desc_res))), desc_res)
    desc_res$Matching[duplicated(desc_res$Matching)] <- ""
  }

  # output
  desc_res
}

#' helper function: sort out a nice report table to summarize binary analysis results
#'
#' @param binobj object from glm()
#' @param weighted_result object res_AB
#' @param eff_measure a string, binary effect measure, could be "OR", "RR", "RD“
#' @param tag a string, by default NULL, if specified, an extra 1st column is created in the output
#'
#' @return a data frame with sample size, incidence rate, estimate of binary effect measure with
#' 95% CI and Wald test of hazard ratio
#'
#' @export

report_table_binary <- function(binobj, weighted_result = NULL, eff_measure = c("OR", "RD", "RR"), tag = NULL) {
  weighted <- ifelse(is.null(weighted_result), FALSE, TRUE)
  # descriptive part
  ARM <- levels(binobj$data$ARM)
  if (!weighted) {
    N <- tapply(binobj$data$USUBJID, binobj$data$ARM, length)
    N.EVNT <- tapply(binobj$data$RESPONSE, binobj$data$ARM, sum)
  } else {
    N <- tapply(binobj$data$weights, binobj$data$ARM, length)
    N.EVNT <- tapply(binobj$data$weights * binobj$data$RESPONSE, binobj$data$ARM, sum)
  }
  N.EVNT.PERC <- round(N.EVNT * 100 / N, 1)
  N.EVNT <- round(N.EVNT, 1)

  # inferential part
  if(!weighted){
    bin_res_est <- coef(binobj)[2]
    bin_res_ci <- confint(binobj,parm = 2, test= "LRT", trace = FALSE)
    bin_res <- c(bin_res_est,bin_res_ci)
    if(eff_measure!="RD") bin_res <- exp(bin_res)
    bin_res <- round(bin_res,2) |> format(nsmall=2)
    bin_res <- paste0(bin_res[1], "[", bin_res[2], ";", bin_res[3], "]")
    bin_pval <- round(summary(binobj)$coefficients[2,4],3) |> format(nsmall=3)
  }else{
    bin_res <- round(unlist(weighted_result[c("est", "ci_l", "ci_u")]), 2) |> format(nsmall = 2)
    bin_res <- paste0(bin_res[1], "[", bin_res[2], ";", bin_res[3], "]")
    bin_pval <- round(weighted_result$pval, 3) |> format(nsmall=3)
  }
  if (bin_pval == "0.000") bin_pval <- "<0.001"

  # assemble the table
  desc_res <- data.frame(
    "treatment" = ARM,
    "N" = N,
    "n.events(%)" = paste0(N.EVNT, "(", format(N.EVNT.PERC, nsmall = 1), ")")
  )

  desc_res <- cbind(desc_res[c(2, 1), ], "xx[95% CI]" = c(bin_res, ""), "p-Value" = c(bin_pval, ""))
  names(desc_res)[ncol(desc_res) - 1] <- paste0(eff_measure, "[95% CI]")
  names(desc_res)[3] <- "n.events(%)" # cbind changes column name

  # add first tag column if applicable
  if (!is.null(tag)) {
    desc_res <- cbind(data.frame(Matching = rep(tag, nrow(desc_res))), desc_res)
    desc_res$Matching[duplicated(desc_res$Matching)] <- ""
  }

  # output
  desc_res
}
