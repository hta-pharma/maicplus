#' Create pseudo IPD given aggregated binary data
#'
#' @param binary_agd a data.frame that take different formats depending on \code{format}
#' @param format a string, "stacked" or "unstacked"
#'
#' @return a data.frame of pseudo binary IPD, with columns USUBJID, ARM, RESPONSE
#' @example inst/examples/get_pseudo_ipd_binary_ex.R
#' @export

get_pseudo_ipd_binary <- function(binary_agd, format = c("stacked", "unstacked")) {
  # pre check
  if (format == "stacked") {
    if (!is.data.frame(binary_agd)) {
      stop("stacked binary_agd should be data.frame with columns 'ARM', 'RESPONSE', 'COUNT'")
    }
    names(binary_agd) <- toupper(names(binary_agd))
    if (!all(c("ARM", "RESPONSE", "COUNT") %in% names(binary_agd))) {
      stop("stacked binary_agd should be data.frame with columns 'ARM', 'Response', 'Count'")
    }
    if (!is.logical(binary_agd$RESPONSE) && !all(toupper(binary_agd$RESPONSE) %in% c("YES", "NO"))) {
      stop("'RESPONSE' column in stacked binary_agd should be either logical vector or character vector of 'Yes'/'No'")
    }
    if (nrow(binary_agd) %% 2 != 0) {
      stop("nrow(binary_agd) is not even number, you may miss to provide 1 level of binary response to certain arm")
    }
  } else if (format == "unstacked") {
    if (!(is.data.frame(binary_agd) || is.matrix(binary_agd))) {
      stop("unstacked binary_agd should be either a 1x2 or 2x2 data frame or matrix")
    }
    if (ncol(binary_agd) != 2 || !nrow(binary_agd) %in% c(1, 2)) {
      stop("unstacked binary_agd should be either a 1x2 or 2x2 data frame or matrix")
    }
    bin_res <- toupper(colnames(binary_agd))
    bin_res <- sort(bin_res)
    if (!(identical(bin_res, c("FALSE", "TRUE")) || identical(bin_res, c("NO", "YES")))) {
      stop("column names of unstacked binary_agd should be either TRUE/FALSE or Yes/No")
    }
  }

  # pre process binary_agd, depending on format
  use_binary_agd <- switch(format,
    "stacked" = {
      names(binary_agd) <- toupper(names(binary_agd))
      if (!is.logical(binary_agd$RESPONSE)) {
        binary_agd$RESPONSE <- toupper(binary_agd$RESPONSE)
        binary_agd$RESPONSE <- binary_agd$RESPONSE == "YES"
      }
      binary_agd
    },
    "unstacked" = {
      trt_names <- rownames(binary_agd)
      bin_res <- toupper(colnames(binary_agd))
      if ("YES" %in% bin_res) {
        bin_res <- ifelse(bin_res == "YES", "TRUE", "FALSE")
        colnames(binary_agd) <- bin_res
      }
      tmpout <- utils::stack(binary_agd)
      tmpout <- cbind(ARM = rep(trt_names, each = 2), tmpout)
      names(tmpout) <- c("ARM", "COUNT", "RESPONSE")
      rownames(tmpout) <- NULL
      tmpout$RESPONSE <- as.logical(tmpout$RESPONSE)
      tmpout
    }
  )

  # create pseudo binary IPD
  use_binary_agd$ARM <- factor(use_binary_agd$ARM, levels = unique(use_binary_agd$ARM))
  n_per_arm <- tapply(use_binary_agd$COUNT, use_binary_agd$ARM, sum)
  n_yes_per_arm <- use_binary_agd$COUNT[use_binary_agd$RESPONSE] # use_binary_agd is already ordered as per factor ARM

  tmpipd <- data.frame(
    USUBJID = NA,
    ARM = unlist(
      mapply(rep, x = levels(use_binary_agd$ARM), each = n_per_arm, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    ),
    RESPONSE = unlist(
      lapply(seq_along(n_per_arm), function(ii) {
        c(rep(TRUE, n_yes_per_arm[ii]), rep(FALSE, n_per_arm[ii] - n_yes_per_arm[ii]))
      })
    )
  )
  tmpipd$USUBJID <- paste0("pseudo_binary_subj_", seq_len(nrow(tmpipd)))

  # output
  tmpipd
}


#' Helper function to summarize outputs from glm fit
#'
#' @param binobj returned object from \code{stats::glm}
#' @param legend label to indicate the binary fit
#' @param weighted logical flag indicating whether weights have been applied in the glm fit
#' @returns A `data.frame` containing a summary of the number of events and subjects in a logistic
#' regression model.
#' @examples
#' data(adrs_sat)
#' pseudo_adrs <- get_pseudo_ipd_binary(
#'   binary_agd = data.frame(
#'     ARM = rep("B", 2),
#'     RESPONSE = c("YES", "NO"),
#'     COUNT = c(280, 120)
#'   ),
#'   format = "stacked"
#' )
#' pseudo_adrs$RESPONSE <- as.numeric(pseudo_adrs$RESPONSE)
#' combined_data <- rbind(adrs_sat[, c("USUBJID", "ARM", "RESPONSE")], pseudo_adrs)
#' combined_data$ARM <- as.factor(combined_data$ARM)
#' binobj_dat <- stats::glm(RESPONSE ~ ARM, combined_data, family = binomial("logit"))
#' glm_makeup(binobj_dat)
#' @export
glm_makeup <- function(binobj, legend = "before matching", weighted = FALSE) {
  arm <- levels(binobj$data$ARM)
  if (!weighted) {
    n <- tapply(binobj$data$USUBJID, binobj$data$ARM, length)
    n_event <- tapply(binobj$data$RESPONSE, binobj$data$ARM, sum)
  } else {
    n <- tapply(binobj$data$weights, binobj$data$ARM, length)
    n_event <- tapply(binobj$data$weights * binobj$data$RESPONSE, binobj$data$ARM, sum)
  }
  data.frame(
    treatment = arm,
    type = legend,
    n = n,
    events = n_event,
    events_pct = n_event * 100 / n
  )
}
