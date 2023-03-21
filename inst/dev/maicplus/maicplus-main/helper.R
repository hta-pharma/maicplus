#' Calculate pooled arm statistics in AgD based on arm-specific statistics (internal use)
#'
#' @param input
#'
#' @return
#' @export
#'
#' @examples
complete_agd <- function(input) {
  input <- as.data.frame(input)
  id <- with(input, which(tolower(treatment) == "total"))
  if (id != nrow(input)) stop("treatment='Total' should be your last row")
  cid <- which(colnames(input) == "pts_n") + 1
  wid <- is.na(input[id, ])
  nid <- grepl("_n$", names(input))
  if (any(wid & nid)) input[id, wid & nid] <- apply(input[1:(id - 1), wid & nid, drop = F], 2, sum)
  mid <- grepl("_mean$", names(input))
  if (any(wid & mid)) {
    ww <- input$pts_n[1:(id - 1)] / sum(input$pts_n[1:(id - 1)])
    input[id, wid & mid] <- ww %*% input[1:(id - 1), wid & mid]
  }
  input
}


#' helper function: transform TTE ADaM data to suitable input for survival R pkg
#'
#' @param dd data frame, ADTTE read via haven::read_sas
#' @param time_scale a character string, 'year', 'month', 'week' or 'day', time unit of median survival time
#'
#' @return a data frame can be used as input to survival::Surv
ext_tte_transfer <- function(dd, time_scale = "month", trt = NULL) {
  timeUnit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)

  if (!time_scale %in% names(timeUnit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")

  if ("CENSOR" %in% names(dd)) {
    dd <- dd[!is.na(dd$CENSOR), ]
    dd$status <- as.numeric(dd$CENSOR == 0)
  }
  if ("EVENT" %in% names(dd)) {
    dd$status <- as.numeric(as.character(dd$EVENT))
  }
  if ("TIME" %in% names(dd)) {
    dd$AVAL <- as.numeric(as.character(dd$TIME))
  }

  dd$time <- dd$AVAL * timeUnit[[time_scale]]
  if (!is.null(trt)) dd$treatment <- trt
  as.data.frame(dd)
}




#' Calculate weights given internal data convention
#'
#' @param ipd data frame, typically ADSL, but include only the individuals that are targeted for the matching
#' @param agd data frame, only one row, either from the arm (arm-specific matching) or from the study (arm-pooled matching)
#'
#' @return a list, similar to \code{\link{cal_weigths}} with an extra element of 'centered_data'
#' @export
#'
msd_matching <- function(ipd, agd) {
  if (!is.data.frame(ipd)) stop("ipd should be a data frame")
  if (!is.data.frame(agd)) stop("agd should be a data frame")
  if (nrow(agd) > 1) stop("only provide AgD for the interested arm or the study (1 row)")

  names(agd) <- toupper(names(agd)) # in case there are sloppy typo with small letters
  em_names <- setdiff(names(agd), c("TREATMENT", "PTS_N"))

  # swap median with proper data in IPD and AgD
  med_id <- grep(pattern = "_MEDIAN$", em_names)
  for (i in seq_along(med_id)) {
    tmp_var <- gsub("_MEDIAN$", "", em_names[i])
    ipd[[em_names[i]]] <- as.numeric(ipd[[tmp_var]] >= agd[[em_names[i]]])
    agd[[em_names[i]]] <- 0.5
  }

  # calculation
  use_em <- ipd[, em_names]
  ext_stat <- unlist(agd[, !names(agd) %in% c("TREATMENT", "PTS_N")])

  use_em <- as.matrix(as.data.frame(use_em))
  use_em <- sweep(use_em, 2, STATS = ext_stat, FUN = "-")

  tmp <- maicplus::cal_weights(EM = use_em)

  # output
  c(tmp, list(centered_data = use_em))
}
